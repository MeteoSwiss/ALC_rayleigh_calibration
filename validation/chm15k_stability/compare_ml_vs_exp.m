%% compare_ml_vs_exp.m
% Report 19: which predicts the CHM15k Rayleigh constant better, the housekeeping ML model of
% estimate_calib_from_housekeeping_v4.m (bagged regression tree on p_calc, temp_int, base,
% laser_pulses, state_optics, cs) or the simple per-unit model C = C_ref * exp(kappa * (T_int - T_ref))?
%
% Data: the 11-year daily archive A:\CHM15k\<stn>\chm_<stn>_1d_01-Jan-2013-31-Dec-2024.mat, cleaned exactly as
% in estimate_calib_from_housekeeping_v4.m. Train = 2013-2023, test = 2024 (the split v4 defines).
%
% Scenario A - housekeeping only (no recent calibration):
%   ML_A  : bagged tree retrained on the training years (same predictors, 30 bagged trees, min leaf 8,
%           i.e. the Regression Learner "Bagged Trees" preset used for the saved model).
%   EXP_A : per optical module (serlom): ln C = a + kappa (T_int - T_ref), kappa from 31-d anomalies and
%           a from the training years of that module (median of the temperature-normalised constant).
%   SAVED : the saved CL_model_bagged_tree_202502 model as-is (it saw 2024 -> in-sample, reference only).
% Scenario B - operational (the last 30 days of calibrations are known):
%   RM_B  : running median of the previous 30 days of Rayleigh constants (Kalman-like baseline).
%   EXP_B : running median of the previous 30 days of exp(-kappa dT)-normalised constants x exp(kappa dT tonight).
%   ML_B  : ML_A x running median of the previous 30 days of C / ML_A (ML corrected for its recent bias).
% Metric on the 2024 nights: robust spread (1.4826 MAD) and median of ln(C_obs / C_pred), in %.

clear;
clc;
close all;
outdir = 'C:\DATA\Projects\202610_CHM15k_stability';
figdir = 'C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration\doc\reports\figs_chm15k_stability';
stations = {'pay','kse','eri','ber','dav','cbw','gal','ino','joy','lei','lin','mhd','sir'};
predictors = {'p_calc','temp_int','base','laser_pulses','state_optics','cs'};
splitDate = datetime(2024, 1, 1);
window = days(30);

%% Load and clean (as estimate_calib_from_housekeeping_v4.m)
T = table();
for s = 1:numel(stations)
    f = fullfile('A:\CHM15k', stations{s}, ['chm_' stations{s} '_1d_01-Jan-2013-31-Dec-2024.mat']);
    S = load(f, 'T');
    Ts = S.T;
    if istimetable(Ts)
        Ts = timetable2table(Ts);
    end
    % Outliers per optical module (serlom), NOT over the whole station record: a module swap is a
    % level step, and a station-wide screen can flag an entire module era as outliers.
    mods_s = unique(Ts.serlom(isfinite(Ts.serlom)));
    for k = 1:numel(mods_s)
        m = Ts.serlom == mods_s(k);
        Ts.Lidar_constant(m) = filloutliers(Ts.Lidar_constant(m), NaN);
    end
    Ts.station = repmat(string(stations{s}), height(Ts), 1);
    keep = [{'Time','Lidar_constant','serlom','station'}, predictors];
    T = [T; Ts(:, keep)]; %#ok<AGROW>
end
T(T.temp_int < 250, :) = [];
T(T.p_calc <= 0, :) = [];
T(T.laser_pulses <= 0, :) = [];
T(T.base > 1e10, :) = [];
T(isnan(T.Lidar_constant), :) = [];
T.Y = T.Lidar_constant / 1e11;
T.Tc = T.temp_int - 273.15;
T.isTrain = T.Time < splitDate;
fprintf('rows: %d train, %d test\n', sum(T.isTrain), sum(~T.isTrain));

%% Scenario A: ML retrained on 2013-2023
tmpl = templateTree('MinLeafSize', 8);
Xtr = T(T.isTrain, predictors);
% Bagging is random: average the prediction of 10 independently seeded ensembles (fixed seeds ->
% reproducible), which also gives the ML its best (lowest-variance) version.
nSeeds = 10;
P = NaN(height(T), nSeeds);
for sd = 1:nSeeds
    rng(sd);
    mdl = fitrensemble(Xtr, T.Y(T.isTrain), 'Method', 'Bag', 'NumLearningCycles', 30, 'Learners', tmpl);
    P(:, sd) = predict(mdl, T(:, predictors));
end
T.ML_A = mean(P, 2);
S = load('C:\Users\hervo\OneDrive\Documents\MATLAB\ALC\CL_model_bagged_tree_202502.mat');
T.SAVED = S.trainedModel.predictFcn(T);

%% Scenario A: per-module exponential model
T.EXP_A = NaN(height(T), 1);
T.kappa = NaN(height(T), 1);
mods = unique(T.serlom);
for k = 1:numel(mods)
    m = T.serlom == mods(k);
    tr = m & T.isTrain;
    if sum(tr) < 40
        continue;
    end
    t = datenum(T.Time(tr));
    y = log(T.Y(tr));
    x = T.Tc(tr);
    ya = y - movmedian_time(t, y, 31);
    xa = x - movmedian_time(t, x, 31);
    ok = isfinite(ya) & isfinite(xa);
    b = polyfit(xa(ok), ya(ok), 1);
    kappa = b(1);
    Tref = median(x, 'omitnan');
    % level from the last year of training of this module (laser ageing makes older years stale)
    last = t >= max(t) - 365;
    a = median(y(last) - kappa * (x(last) - Tref), 'omitnan');
    T.EXP_A(m) = exp(a + kappa * (T.Tc(m) - Tref));
    T.kappa(m) = kappa;
end

%% Scenario B: operational predictions from the previous 30 days
T.RM_B = NaN(height(T), 1);
T.EXP_B = NaN(height(T), 1);
T.ML_B = NaN(height(T), 1);
for k = 1:numel(mods)
    m = find(T.serlom == mods(k));
    if numel(m) < 40
        continue;
    end
    tt = T.Time(m);
    y = T.Y(m);
    kappa = T.kappa(m(1));
    if ~isfinite(kappa)
        kappa = 0;
    end
    yn = y .* exp(-kappa * T.Tc(m));
    ratio = y ./ T.ML_A(m);
    for i = 1:numel(m)
        prev = tt < tt(i) & tt >= tt(i) - window;
        if sum(prev) < 5
            continue;
        end
        T.RM_B(m(i)) = median(y(prev), 'omitnan');
        T.EXP_B(m(i)) = median(yn(prev), 'omitnan') * exp(kappa * T.Tc(m(i)));
        T.ML_B(m(i)) = T.ML_A(m(i)) * median(ratio(prev), 'omitnan');
    end
end

%% Scores on the 2024 test nights
models = {'ML_A','EXP_A','SAVED','RM_B','EXP_B','ML_B'};
labels = {'ML (HK only, retrained)','exp(\kappa T) (HK only)','ML saved 202502 (in-sample!)', ...
          'running median 30 d','running median x exp(\kappa \DeltaT)','ML x 30-d bias'};
te = ~T.isTrain;
res = table();
for j = 1:numel(models)
    p = T.(models{j});
    ok = te & isfinite(p) & p > 0 & isfinite(T.Y);
    r = 100 * log(T.Y(ok) ./ p(ok));
    row = table(string(models{j}), string(labels{j}), sum(ok), median(r), 1.4826 * median(abs(r - median(r))), ...
        sqrt(mean(r .^ 2)), 'VariableNames', {'model','label','n','median_pct','robust_sd_pct','rmse_pct'});
    res = [res; row]; %#ok<AGROW>
end
% same nights for every model (fair comparison)
okAll = te & isfinite(T.Y);
for j = 1:numel(models)
    okAll = okAll & isfinite(T.(models{j})) & T.(models{j}) > 0;
end
res.robust_sd_common = NaN(height(res), 1);
res.n_common = repmat(sum(okAll), height(res), 1);
for j = 1:numel(models)
    r = 100 * log(T.Y(okAll) ./ T.(models{j})(okAll));
    res.robust_sd_common(j) = 1.4826 * median(abs(r - median(r)));
end
disp(res);
writetable(res, fullfile(outdir, 'ml_vs_exp_scores.csv'));

%% Per station (common nights)
st = unique(T.station(okAll));
per = table();
for k = 1:numel(st)
    m = okAll & T.station == st(k);
    row = table(st(k), sum(m), 'VariableNames', {'station','n'});
    for j = 1:numel(models)
        r = 100 * log(T.Y(m) ./ T.(models{j})(m));
        row.(models{j}) = 1.4826 * median(abs(r - median(r)));
    end
    per = [per; row]; %#ok<AGROW>
end
disp(per);
writetable(per, fullfile(outdir, 'ml_vs_exp_per_station.csv'));
writetable(T(:, [{'Time','station','serlom','Y','Tc','kappa'}, models]), fullfile(outdir, 'ml_vs_exp_predictions.csv'));

%% Figure (landscape)
fig = figure('Position', [50 50 1900 560], 'Color', 'w');
tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
b = bar(categorical(labels, labels), res.robust_sd_common);
b.FaceColor = [0.3 0.5 0.8];
ylabel('robust spread of ln(C_{obs}/C_{pred}) on 2024 nights [%]');
title(sprintf('Prediction error, same %d test nights', res.n_common(1)));
grid on;
nexttile;
cols = lines(numel(models));
xs = 1:height(per);
for j = [1 2 5 6]
    plot(xs, per.(models{j}), 'o-', 'Color', cols(j, :), 'LineWidth', 1.4, 'DisplayName', labels{j});
    hold on;
end
xticks(xs);
xticklabels(cellstr(per.station));
ylabel('robust spread [%]');
title('Per station (2024 test nights)');
legend('Location', 'best');
grid on;
nexttile;
m = okAll & T.station == "pay";
plot(T.Time(m), T.Y(m), 'k.', 'MarkerSize', 8, 'DisplayName', 'Rayleigh C (observed)');
hold on;
plot(T.Time(m), T.ML_A(m), '-', 'Color', cols(1, :), 'DisplayName', labels{1});
plot(T.Time(m), T.EXP_A(m), '-', 'Color', cols(2, :), 'DisplayName', labels{2});
plot(T.Time(m), T.EXP_B(m), '-', 'Color', cols(5, :), 'LineWidth', 1.4, 'DisplayName', labels{5});
ylabel('C_L / 1e11');
title('Payerne 2024 (test year)');
legend('Location', 'best');
grid on;
exportgraphics(fig, fullfile(figdir, 'q8_ml_vs_exp.png'), 'Resolution', 150);
fprintf('ML_VS_EXP_DONE\n');

%% Local functions
function out = movmedian_time(t, x, win)
% centred running median over +-win/2 days on irregular sampling (t in datenum)
out = NaN(size(x));
for i = 1:numel(t)
    m = abs(t - t(i)) <= win / 2 & isfinite(x);
    if sum(m) >= 5
        out(i) = median(x(m));
    end
end
end
