% ===== Provenance (paper code archive) =====
% Paper Figure A5-1: dark (covered-telescope) measurement of 12 May 2026 for the Payerne CL31, CL61 and CHM15k
% ("Figure 2" section, lines ~299-360 of the original file).
% CAVEAT (2026-10-20): this saved version does NOT reproduce the paper figure exactly: the figure title says
% 14:05 and its y-axis stops at 8000 m, whereas this file has t2 = 14:50 and range_top = 15000. The figure
% (in the manuscript since v0.5, 2026-06-23) came from an earlier, unsaved edit. Conventions also differ from
% Figure A5-2: CHM15k divided by a hard-coded 5e11 (operational constant that day 4.73e11) and CL31 by a
% nominal 1e8 (operational cloud constant 2.83e7).
% Source: C:/Users/hervo/OneDrive/Documents/MATLAB/ALC (not under git), last modified 2026-07-03.
% Copied here unchanged on 2026-10-20.
% ==========================================
clearvars;close all;clc;

%% input
t1 = datetime(2026,5,12,9,35,0);
t2 = datetime(2026,5,12,14,50,0);

% t1 = datetime(2026,5,26,11,45,0);
% t2 = datetime(2026,5,27,13,15,0);
% 
% t1 = datetime(2026,6,9,9,20,0);
% t2 = datetime(2026,6,9,11,50,0);


% t1 = datetime(2026,5,12,1,35,0);
% t2 = datetime(2026,5,12,14,50,0);

% CHM15k calibration constant for the Payerne unit.
% beta_att [m^-1 sr^-1] = beta_raw / CHM_CAL
CHM_CAL = 5e11;% Update this if the unit is re-calibrated or replaced.

%% get data
cl31 = get_cl31_from_files( 'pay', t1, t2);
cl61 = get_cl61_from_files( 'pay', t1, t2);
chm  = get_chm15k_from_files('pay', t1, t2);


%% background
figure
subplot(3,1,1)
plot(cl31.dtime,cl31.bckgrd_rcs_910)
grid on
ylabel('Background RCS')
title('CL31: bckgrd\_rcs\_910')

subplot(3,1,2)
plot(chm.dtime,chm.base)
grid on
ylabel('Baseline')
title('CHM15k: base')

subplot(3,1,3)
plot(cl61.dtime,cl61.beta_att_noise_level)
grid on
ylabel('Noise level [unitless]')
xlabel('Time (UTC)')
title('CL61: beta\_att\_noise\_level')
%% Bring everything to a common scale: attenuated backscatter in Mm^-1 sr^-1
% References:
%   CL31  : Vaisala User's Guide M210482EN; E-PROFILE / Raw2L1 convention.
%           rcs_910 in the bulletin is the range-corrected attenuated
%           backscatter stored as integers in "100000 km^-1 sr^-1" units
%           (i.e. 1 LSB = 1e-5 km^-1 sr^-1 = 1e-8 m^-1 sr^-1).
%           Conversion: beta_att [m^-1 sr^-1] = rcs_910 * 1e-8.
%   CL61  : NetCDF beta_att and beta_att_noise_level are already in
%           m^-1 sr^-1 (CF-compliant).
%   CHM15k: NetCDF beta_raw is range-corrected and calibration-normalised.
%           Dividing by the calibration constant (variable name depends
%           on the file: c_cal, calibration_factor, scaling, ...) gives
%           beta_att in m^-1 sr^-1.

UNIT_SCALE = 1e6;   % m^-1 sr^-1 -> Mm^-1 sr^-1 for display

% CL31
r_cl31 = cl31.range;
beta_cl31 = cl31.rcs_910 * 1e-8 * UNIT_SCALE;   % -> Mm^-1 sr^-1

% CL61
r_cl61 = cl61.range;
beta_cl61 = cl61.beta_att * UNIT_SCALE;          % -> Mm^-1 sr^-1

% CHM15k - apply the user-defined calibration constant.
% beta_raw is range-corrected and calibration-normalised. Dividing by the
% per-instrument calibration constant gives beta_att in m^-1 sr^-1.
% For the Payerne CHM15k unit, CHM_CAL = 5e11 (set at the top of the script).
% If the NetCDF file ships a calibration variable, we report it for
% comparison but keep using CHM_CAL as the authoritative value.
r_chm = chm.range;
if size(chm.beta_raw, 1) == numel(r_chm)
    beta_chm_rc = chm.beta_raw';
else
    beta_chm_rc = chm.beta_raw;
end

chm_cal_candidates = {'c_cal', 'calibration_factor', 'scaling', ...
                      'calibration_const'};
for i = 1:numel(chm_cal_candidates)
    nm = chm_cal_candidates{i};
    if isfield(chm, nm) && ~isempty(chm.(nm))
        fprintf('CHM15k file ships %s = %.3e (using CHM_CAL = %.3e from script).\n', ...
                nm, double(chm.(nm)(1)), CHM_CAL);
        break;
    end
end
fprintf('Using CHM15k calibration constant: CHM_CAL = %.3e (Payerne unit)\n', CHM_CAL);
beta_chm = (beta_chm_rc / CHM_CAL) * UNIT_SCALE;  % -> Mm^-1 sr^-1

%% Summary stats per instrument
med_cl31 = median(beta_cl31, 1, 'omitnan');
p25_cl31 = prctile(beta_cl31, 25, 1);
p75_cl31 = prctile(beta_cl31, 75, 1);
std_cl31 = std(beta_cl31, 0, 1, 'omitnan');
rms_cl31 = sqrt(mean(beta_cl31.^2, 1, 'omitnan'));

med_cl61 = median(beta_cl61, 1, 'omitnan');
p25_cl61 = prctile(beta_cl61, 25, 1);
p75_cl61 = prctile(beta_cl61, 75, 1);
std_cl61 = std(beta_cl61, 0, 1, 'omitnan');
rms_cl61 = sqrt(mean(beta_cl61.^2, 1, 'omitnan'));

med_chm  = median(beta_chm, 1, 'omitnan');
p25_chm  = prctile(beta_chm, 25, 1);
p75_chm  = prctile(beta_chm, 75, 1);
std_chm  = std(beta_chm, 0, 1, 'omitnan');
rms_chm  = sqrt(mean(beta_chm.^2, 1, 'omitnan'));

%% Colors and labels
col_cl31 = [0.00 0.45 0.74];
col_cl61 = [0.85 0.33 0.10];
col_chm  = [0.47 0.67 0.19];

tstr = sprintf('Dark measurement comparison - Payerne %s to %s', ...
               datestr(t1,'yyyy-mm-dd HH:MM'), datestr(t2,'HH:MM UTC'));

clim_dark = [-0.5 0.5];
range_top = 15000;

%% Figure 1 : time-height panels
bigfig;

subplot(3,1,1)
pcolor(cl31.dtime, r_cl31, beta_cl31'); shading flat
cb = colorbar; cb.Label.String = '\beta_{att} [Mm^{-1} sr^{-1}]';
clim(clim_dark)
ylabel('Range [m]'); title('CL31  (910 nm)')
ylim([0 range_top])
grid on; box on

subplot(3,1,2)
pcolor(cl61.dtime, r_cl61, beta_cl61'); shading flat
cb = colorbar; cb.Label.String = '\beta_{att} [Mm^{-1} sr^{-1}]';
clim(clim_dark)
ylabel('Range [m]'); title('CL61  (910 nm)')
ylim([0 range_top])
grid on; box on

subplot(3,1,3)
pcolor(chm.dtime, r_chm, beta_chm'); shading flat
cb = colorbar; cb.Label.String = '\beta_{att} [Mm^{-1} sr^{-1}]';
clim(clim_dark)
ylabel('Range [m]'); xlabel('Time (UTC)'); title('CHM15k  (1064 nm)')
ylim([0 range_top])
grid on; box on

sgtitle(tstr, 'FontSize', 13, 'FontWeight', 'bold')

%% Plot CL61 with Median
% Step 1: median filter to reject single-gate/single-shot spikes (robust, no smearing of real trend)
beta_medfilt = medfilt2(beta_cl61'./r_cl61./r_cl61, [9 5]);   % [range_gates, time_profiles]

% Step 2: 2-D smoothing to reveal the slowly varying offset
beta_smooth = smoothdata2(beta_medfilt, 'gaussian', {21 61});  % [range_window, time_window]

figure
pcolor(cl61.dtime, r_cl61, beta_medfilt);
shading flat
cb = colorbar; cb.Label.String = '\beta_{att} [Mm^{-1} sr^{-1}]';
clim([-0.1 0.1] *1e-8)
ylabel('Range [m]'); title('CL61  (910 nm) - offset (filtered)')
ylim([0 range_top])
grid on; box on
%%
beta_medfilt = medfilt2(beta_cl61'./r_cl61./r_cl61, [9 5]);   % [range_gates, time_profiles]

    figure
plot(smoothdata(median(beta_medfilt,2),"rlowess",50), r_cl61)
xlim([-0.2 0.2]*1e-8)
grid on
hold on

% Fit beta (x) as a linear function of range (y), for range > 2000 m
r_fit = {[0,2000],[2000 12000],[12000 15000]};
for i = 1:3
beta_profile = smoothdata(median(beta_medfilt(r_cl61>r_fit{i}(1) & r_cl61 <r_fit{i}(2),:),2),"rlowess",50);

range_profile = r_cl61(r_cl61>r_fit{i}(1) & r_cl61 <r_fit{i}(2));

f = polyfit(range_profile, beta_profile, 1);   % beta ≈ f(1)*range + f(2)

beta_fit = polyval(f, r_fit{i});

plot(beta_fit, r_fit{i}, '-', 'LineWidth', 1.5)

end


legend('Median profile (smoothed)', sprintf('Linear fit: %.3g\\cdotr + %.3g', f(1), f(2)), 'Location','best')


%% Plot CHM with Median
% Step 1: median filter to reject single-gate/single-shot spikes (robust, no smearing of real trend)
beta_medfilt = medfilt2(beta_chm'./r_chm./r_chm, [9 5])*(3e11);   % [range_gates, time_profiles]


figure
pcolor(chm.dtime, r_chm, beta_medfilt);
shading flat
cb = colorbar; cb.Label.String = '\beta_{att} [Mm^{-1} sr^{-1}]';
clim([-0.1 0.1] *1e-8)
ylabel('Range [m]'); title('CHM  (1064 nm) -P')
ylim([0 range_top])
grid on; box on

%% Histogram of median(beta/r^2) per 1 km altitude layer

% r_chm assumed to be a range vector (m); make it a column so it
% broadcasts correctly against beta_chm' (range x time)
r_chm = r_chm(:);

edges   = (-10:0.25:10)*1e-9;
z_edges = 0:2500:ceil(max(r_chm)/1000)*1000;
figure
hold on
for k = 1:numel(z_edges)-1
    z_lo = z_edges(k);
    z_hi = z_edges(k+1);
    idx  = r_chm >= z_lo & r_chm < z_hi;

    if ~any(idx)
        continue
    end

    beta_norm = beta_chm(:,idx)'./r_chm(idx)./r_chm(idx);
    med_beta  = median(beta_norm, 2, 'omitnan');

    
    histogram(med_beta, edges,'displayname',sprintf('%d - %d m', z_lo, z_hi)')
    xlabel('median(\beta/r^2)  [m^{-1} sr^{-1}]')
    ylabel('Count')
end
%%

 
figure
subplot(3,2,[1 3])
tmpx = beta_chm'./r_chm./r_chm;
tmpy = repmat(r_chm,1,length(chm.dtime));
binscatter(tmpx(:),tmpy(:))
hold on
plot(median(beta_chm'./r_chm./r_chm,2,'omitnan'), r_chm,'.k')
xlim([-2.5 2.5]*1e-7)
title('chm15k')
ylim([0 15000])
xlabel('\beta/r^2  [m^{-1} sr^{-1}]')
ylabel('Range [m]')
grid on

subplot(3,2,5)
histogram(median(beta_chm'./r_chm./r_chm,2,'omitnan'), (-2.5:0.1:2.5)*1e-7)
xlabel('median(\beta/r^2)  [m^{-1} sr^{-1}]')
ylabel('Count')
grid on
c = colorbar;

subplot(3,2,[2 4])
tmpx = beta_cl61'./r_cl61./r_cl61;
tmpy = repmat(r_cl61,1,length(cl61.dtime));
binscatter(tmpx(:),tmpy(:))
hold on
plot(median(beta_cl61'./r_cl61./r_cl61,2,'omitnan'), r_cl61,'.k')
xlim([-2.5 2.5]*1e-7)
title('CL61')
ylim([0 15000])
xlabel('\beta/r^2  [m^{-1} sr^{-1}]')
ylabel('Range [m]')
grid on

subplot(3,2,6)
histogram(median(beta_cl61'./r_cl61./r_cl61,2,'omitnan'), (-2.5:0.1:2.5)*1e-7)
xlabel('median(\beta/r^2)  [m^{-1} sr^{-1}]')
ylabel('Count')
grid on
c = colorbar;



sgtitle(sprintf('Terminal hood dark measurements \n  %s to %s ', ...
    datestr(chm.dtime(1),'dd-mmm-yyyy HH:MM'), datestr(cl61.dtime(end),'dd-mmm-yyyy HH:MM')))
%%
% Fit beta (x) as a linear function of range (y), for range > 2000 m
r_fit = {[0,2000],[2000 12000],[12000 15000]};
for i = 1:3
    beta_profile = smoothdata(median(beta_medfilt(r_chm>r_fit{i}(1) & r_chm <r_fit{i}(2),:),2),"rlowess",50);

    range_profile = r_chm(r_chm>r_fit{i}(1) & r_chm <r_fit{i}(2));

    f = polyfit(range_profile, beta_profile, 1);   % beta ≈ f(1)*range + f(2)

    beta_fit = polyval(f, r_fit{i});

    plot(beta_fit, r_fit{i}, '-', 'LineWidth', 1.5)

end
%% Figure 2 : mean profile (IQR shading) + noise vs range
bigfig;

subplot(1,3,1)
hold on
fill([p25_cl31 fliplr(p75_cl31)], [r_cl31; flipud(r_cl31)]', ...
     col_cl31, 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(med_cl31, r_cl31, '-', 'Color', col_cl31, 'LineWidth', 1.5, 'DisplayName','CL31');

fill([p25_cl61 fliplr(p75_cl61)], [r_cl61; flipud(r_cl61)]', ...
     col_cl61, 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(med_cl61, r_cl61, '-', 'Color', col_cl61, 'LineWidth', 1.5, 'DisplayName','CL61');

fill([p25_chm fliplr(p75_chm)], [r_chm; flipud(r_chm)]', ...
     col_chm, 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(med_chm, r_chm, '-', 'Color', col_chm, 'LineWidth', 1.5, 'DisplayName','CHM15k');

xline(0, 'k:', 'HandleVisibility','off')
hold off
box on; grid on
xlabel('\beta_{att} [Mm^{-1} sr^{-1}]')
ylabel('Range [m]')
title('Median profile (shaded: IQR 25/75 %)')
xlim([-0.5 0.5])
ylim([0 range_top])
legend('Location','best')

subplot(1,3,2)
hold on
plot(std_cl31, r_cl31, '-', 'Color', col_cl31, 'LineWidth', 1.5, 'DisplayName','CL31');
plot(std_cl61, r_cl61, '-', 'Color', col_cl61, 'LineWidth', 1.5, 'DisplayName','CL61');
plot(std_chm,  r_chm,  '-', 'Color', col_chm,  'LineWidth', 1.5, 'DisplayName','CHM15k');
hold off
box on; grid on
set(gca,'xscale','log')
xlabel('\sigma(\beta_{att}) [Mm^{-1} sr^{-1}]')
ylabel('Range [m]')
title('Temporal noise (std)')
xlim([1e-4 1e2])
ylim([0 range_top])
legend('Location','best')

subplot(1,3,3)
hold on
plot(rms_cl31, r_cl31, '-', 'Color', col_cl31, 'LineWidth', 1.5, 'DisplayName','CL31');
plot(rms_cl61, r_cl61, '-', 'Color', col_cl61, 'LineWidth', 1.5, 'DisplayName','CL61');
plot(rms_chm,  r_chm,  '-', 'Color', col_chm,  'LineWidth', 1.5, 'DisplayName','CHM15k');
hold off
box on; grid on
set(gca,'xscale','log')
xlabel('RMS(\beta_{att}) [Mm^{-1} sr^{-1}]')
ylabel('Range [m]')
title('RMS = \surd(bias^2 + \sigma^2)')
xlim([1e-4 1e2])
ylim([0 range_top])
legend('Location','best')

sgtitle(tstr, 'FontSize', 13, 'FontWeight', 'bold')

%% Figure 3 : near-range zoom
figure;
hold on
plot(med_cl31, r_cl31, '-', 'Color', col_cl31, 'LineWidth', 1.5, 'DisplayName','CL31');
plot(med_cl61, r_cl61, '-', 'Color', col_cl61, 'LineWidth', 1.5, 'DisplayName','CL61');
plot(med_chm,  r_chm,  '-', 'Color', col_chm,  'LineWidth', 1.5, 'DisplayName','CHM15k');
xline(0, 'k:', 'HandleVisibility','off')
hold off
box on; grid on
xlabel('\beta_{att} [Mm^{-1} sr^{-1}]')
ylabel('Range [m]')
title(['Near-range median profile - ' tstr])
xlim([-0.3 0.3])
ylim([0 1500])
legend('Location','best')

%% Summary statistics (printed)
band = [500 3000];
fprintf('\n=== Dark-measurement summary (mean over %g-%g m) ===\n', band(1), band(2));

mask_31_print = r_cl31 >= band(1) & r_cl31 <= band(2);
mask_61_print = r_cl61 >= band(1) & r_cl61 <= band(2);
mask_ch_print = r_chm  >= band(1) & r_chm  <= band(2);

fprintf('%-8s  %-12s  %-12s  %-12s  %-7s  %-7s\n', ...
        'Instr', 'bias[Mm-1]', 'std[Mm-1]', 'RMS[Mm-1]', 'Nprof', 'Nrange');
fprintf('%-8s  %-12.4f  %-12.4f  %-12.4f  %-7d  %-7d\n', 'CL31', ...
        mean(med_cl31(mask_31_print), 'omitnan'), ...
        mean(std_cl31(mask_31_print), 'omitnan'), ...
        mean(rms_cl31(mask_31_print), 'omitnan'), ...
        size(beta_cl31,1), nnz(mask_31_print));
fprintf('%-8s  %-12.4f  %-12.4f  %-12.4f  %-7d  %-7d\n', 'CL61', ...
        mean(med_cl61(mask_61_print), 'omitnan'), ...
        mean(std_cl61(mask_61_print), 'omitnan'), ...
        mean(rms_cl61(mask_61_print), 'omitnan'), ...
        size(beta_cl61,1), nnz(mask_61_print));
fprintf('%-8s  %-12.4f  %-12.4f  %-12.4f  %-7d  %-7d\n', 'CHM15k', ...
        mean(med_chm(mask_ch_print), 'omitnan'), ...
        mean(std_chm(mask_ch_print), 'omitnan'), ...
        mean(rms_chm(mask_ch_print), 'omitnan'), ...
        size(beta_chm,1), nnz(mask_ch_print));
fprintf('=====================================================\n\n');

%% Noise diagnostics
% DARK MEASUREMENT: telescopes are covered, so the detector sees no
% atmospheric signal at any range. Every range gate is independent dark
% noise. We use the full range (instead of restricting to high altitude)
% to maximize the number of independent samples per averaging window.

mask_31 = true(size(r_cl31));
mask_61 = true(size(r_cl61));
mask_ch = true(size(r_chm));

noise_per_prof_31 = std(beta_cl31(:, mask_31), 0, 2, 'omitnan');
noise_per_prof_61 = std(beta_cl61(:, mask_61), 0, 2, 'omitnan');
noise_per_prof_ch = std(beta_chm( :, mask_ch), 0, 2, 'omitnan');

%% Pull temperatures
% CL31
T_cl31 = NaN(size(cl31.dtime));
if isfield(cl31, 'laser_temp')
    T_cl31 = cl31.laser_temp(:);
end

% CL61
T_cl61 = NaN(size(cl61.dtime));
cl61_temp_candidates = { ...
    'monitoring_laser_temperature', ...
    'monitoring_transmitter_enclosure_temperature', ...
    'monitoring_internal_temperature'};
T_cl61_name = '';
for i = 1:numel(cl61_temp_candidates)
    if isfield(cl61, cl61_temp_candidates{i})
        T_cl61 = cl61.(cl61_temp_candidates{i})(:);
        T_cl61_name = cl61_temp_candidates{i};
        break;
    end
end

% CHM15k
T_chm  = NaN(size(chm.dtime));
chm_temp_candidates = {'temp_int', 'temp_det', 'temp_lom', 'temp_ext'};
T_chm_name = '';
for i = 1:numel(chm_temp_candidates)
    if isfield(chm, chm_temp_candidates{i})
        T_chm = chm.(chm_temp_candidates{i})(:);
        T_chm_name = chm_temp_candidates{i};
        break;
    end
end

% CHM temp_int is reported in Kelvin in the standard NetCDF; convert to degC
if ~isempty(T_chm_name) && any(isfinite(T_chm)) && median(T_chm, 'omitnan') > 200
    T_chm = T_chm - 273.15;
end

if isempty(T_cl61_name)
    warning('No CL61 temperature variable found.');
else
    fprintf('Using CL61 temperature variable: %s\n', T_cl61_name);
end
if isempty(T_chm_name)
    warning('No CHM15k temperature variable found.');
else
    fprintf('Using CHM15k temperature variable: %s\n', T_chm_name);
end

%% Figure 4 : noise vs temperature
bigfig;

subplot(1,3,1)
plot(T_cl31, noise_per_prof_31, '.', 'MarkerSize', 6);
box on; grid on
xlabel('Laser temperature [\circC]')
ylabel('\sigma(\beta_{att}) full range  [Mm^{-1} sr^{-1}]')
title('CL31: noise vs laser temp')
r = local_corr(T_cl31, noise_per_prof_31);
text(0.02, 0.97, sprintf('r = %+.2f', r), 'Units','normalized', ...
     'VerticalAlignment','top', 'FontWeight','bold');

subplot(1,3,2)
plot(T_cl61, noise_per_prof_61, '.', 'MarkerSize', 6);
box on; grid on
xlabel('Temperature [\circC]')
ylabel('\sigma(\beta_{att}) [Mm^{-1} sr^{-1}]')
if isempty(T_cl61_name)
    title('CL61: NO temperature variable found')
else
    title(sprintf('CL61: noise vs %s', strrep(T_cl61_name,'_','\_')))
end
r = local_corr(T_cl61, noise_per_prof_61);
text(0.02, 0.97, sprintf('r = %+.2f', r), 'Units','normalized', ...
     'VerticalAlignment','top', 'FontWeight','bold');

subplot(1,3,3)
plot(T_chm, noise_per_prof_ch, '.', 'MarkerSize', 6);
box on; grid on
xlabel('Internal temperature [\circC]')
ylabel('\sigma(\beta_{att}) [Mm^{-1} sr^{-1}]')
if isempty(T_chm_name)
    title('CHM15k: NO temperature variable found')
else
    title(sprintf('CHM15k: noise vs %s', strrep(T_chm_name,'_','\_')))
end
r = local_corr(T_chm, noise_per_prof_ch);
text(0.02, 0.97, sprintf('r = %+.2f', r), 'Units','normalized', ...
     'VerticalAlignment','top', 'FontWeight','bold');

sgtitle(['Noise vs internal temperature - ' tstr], 'FontWeight','bold')

%% Figure 5 : time series overlay of noise and temperature
bigfig;

subplot(3,1,1)
yyaxis left
plot(cl31.dtime, noise_per_prof_31, '.', 'Color', col_cl31, 'MarkerSize', 6);
ylabel('\sigma(\beta_{att})')
yyaxis right
plot(cl31.dtime, T_cl31, '-k', 'LineWidth', 1)
ylabel('Laser temp [\circC]')
grid on; box on
title('CL31')

subplot(3,1,2)
yyaxis left
plot(cl61.dtime, noise_per_prof_61, '.', 'Color', col_cl61, 'MarkerSize', 6);
ylabel('\sigma(\beta_{att})')
yyaxis right
plot(cl61.dtime, T_cl61, '-k', 'LineWidth', 1)
ylabel('Temperature [\circC]')
grid on; box on
if isempty(T_cl61_name)
    title('CL61')
else
    title(sprintf('CL61 (%s)', strrep(T_cl61_name,'_','\_')))
end

subplot(3,1,3)
yyaxis left
plot(chm.dtime, noise_per_prof_ch, '.', 'Color', col_chm, 'MarkerSize', 6);
ylabel('\sigma(\beta_{att})')
yyaxis right
plot(chm.dtime, T_chm, '-k', 'LineWidth', 1)
ylabel('Internal temp [\circC]')
grid on; box on
if isempty(T_chm_name)
    title('CHM15k')
else
    title(sprintf('CHM15k (%s)', strrep(T_chm_name,'_','\_')))
end
xlabel('Time (UTC)')

sgtitle(['Noise and temperature time series - ' tstr], 'FontWeight','bold')

%% Figure 6 : CL61 instrument-reported noise vs empirical std
% Vaisala beta_att_noise_level is documented in the CL61 NetCDF as
% "a unitless number describing the noise level of the attenuated
%  volume backscatter coefficient" - it is NOT a noise in m^-1 sr^-1
% despite the variable name. We plot it on its own (unitless) axis and
% the empirical sigma of beta_att over a high-altitude noise band on a
% second axis. They are expected to *track in time* but cannot be
% compared in absolute units.
noise_band = [4000 7000];
mask_61_noise = r_cl61 >= noise_band(1) & r_cl61 <= noise_band(2);
emp_noise_61_band = std(beta_cl61(:, mask_61_noise), 0, 2, 'omitnan');

if isfield(cl61, 'beta_att_noise_level')
    figure;
    yyaxis left
    plot(cl61.dtime, cl61.beta_att_noise_level, '-', ...
         'Color', col_cl61, 'LineWidth', 1.2);
    ylabel('Vaisala beta\_att\_noise\_level  [Unitless]')
    yyaxis right
    plot(cl61.dtime, emp_noise_61_band, '.k', 'MarkerSize', 6);
    ylabel(sprintf('Empirical \\sigma over %g-%g m  [Mm^{-1} sr^{-1}]', ...
                   noise_band(1), noise_band(2)))
    grid on; box on
    xlabel('Time (UTC)')
    title(['CL61: instrument-reported vs empirical noise - ' tstr])
end

%% Allan analysis : noise vs averaging duration
% Two complementary metrics are computed per averaging time tau:
%   (1) Spatial-spread metric (local_allan_noise): per tau-window, the std
%       across range gates of the window time-mean, averaged over windows.
%       This is the original "spread of the tau-averaged dark field" view.
%   (2) Time-domain Allan deviation (local_allan_time): the classic
%       non-overlapping two-sample deviation on a range-collapsed time
%       series. Its minimum marks the optimal averaging time before drift
%       (e.g. the diurnal temperature cycle) starts to dominate over white
%       detector noise. With a >24 h record this knee becomes resolvable.
dt_cl31 = seconds(median(diff(cl31.dtime)));
dt_cl61 = seconds(median(diff(cl61.dtime)));
dt_chm  = seconds(median(diff(chm.dtime)));

fprintf('Sampling rates:  CL31=%.1f s   CL61=%.1f s   CHM=%.1f s\n', ...
        dt_cl31, dt_cl61, dt_chm);

% Span the data: cap tau at ~1/8 of the shortest record so every instrument
% keeps at least ~8 non-overlapping segments at the longest tau.
T_span = min([seconds(cl31.dtime(end)-cl31.dtime(1)), ...
              seconds(cl61.dtime(end)-cl61.dtime(1)), ...
              seconds(chm.dtime(end) -chm.dtime(1))]);
tau_max = T_span / 8;
fprintf('Record span (shortest) = %.1f h  ->  tau_max = %.1f h\n', ...
        T_span/3600, tau_max/3600);

tau_sec = unique(round(logspace(log10(min([dt_cl31 dt_cl61 dt_chm])), ...
                                log10(tau_max), 45)));

% Range-collapsed time series for the time-domain Allan deviation.
ts_cl31 = mean(beta_cl31(:, mask_31), 2, 'omitnan');
ts_cl61 = mean(beta_cl61(:, mask_61), 2, 'omitnan');
ts_chm  = mean(beta_chm( :, mask_ch), 2, 'omitnan');

allan_cl31 = nan(size(tau_sec));
allan_cl61 = nan(size(tau_sec));
allan_chm  = nan(size(tau_sec));
nwin_cl31  = nan(size(tau_sec));
nwin_cl61  = nan(size(tau_sec));
nwin_chm   = nan(size(tau_sec));

adev_cl31  = nan(size(tau_sec));
adev_cl61  = nan(size(tau_sec));
adev_chm   = nan(size(tau_sec));
nseg_cl31  = nan(size(tau_sec));
nseg_cl61  = nan(size(tau_sec));
nseg_chm   = nan(size(tau_sec));

for k = 1:numel(tau_sec)
    tau = tau_sec(k);
    [allan_cl31(k), nwin_cl31(k)] = local_allan_noise(beta_cl31(:, mask_31), dt_cl31, tau);
    [allan_cl61(k), nwin_cl61(k)] = local_allan_noise(beta_cl61(:, mask_61), dt_cl61, tau);
    [allan_chm(k),  nwin_chm(k) ] = local_allan_noise(beta_chm( :, mask_ch), dt_chm,  tau);

    [adev_cl31(k), nseg_cl31(k)] = local_allan_time(ts_cl31, dt_cl31, tau);
    [adev_cl61(k), nseg_cl61(k)] = local_allan_time(ts_cl61, dt_cl61, tau);
    [adev_chm(k),  nseg_chm(k) ] = local_allan_time(ts_chm,  dt_chm,  tau);
end

% Require enough segments for a meaningful estimate at long tau.
min_seg = 8;
allan_cl31(nwin_cl31 < min_seg) = NaN;
allan_cl61(nwin_cl61 < min_seg) = NaN;
allan_chm( nwin_chm  < min_seg) = NaN;
adev_cl31(nseg_cl31 < min_seg) = NaN;
adev_cl61(nseg_cl61 < min_seg) = NaN;
adev_chm( nseg_chm  < min_seg) = NaN;

% Locate the drift knee (minimum of the time-domain Allan deviation).
[tau_opt_31, adev_min_31] = local_allan_min(tau_sec, adev_cl31);
[tau_opt_61, adev_min_61] = local_allan_min(tau_sec, adev_cl61);
[tau_opt_ch, adev_min_ch] = local_allan_min(tau_sec, adev_chm);

%% Figure 7 : Allan analysis vs averaging duration
xt_sec   = [10 30 60 300 600 1800 3600 3*3600 6*3600 12*3600 24*3600];
xt_label = {'10s','30s','1min','5min','10min','30min','1h','3h','6h','12h','24h'};

bigfig;

subplot(1,2,1)
loglog(tau_sec, allan_cl31, '-o', 'Color', col_cl31, 'LineWidth',1.5, ...
       'MarkerFaceColor', col_cl31, 'MarkerSize', 5, 'DisplayName','CL31');
hold on
loglog(tau_sec, allan_cl61, '-s', 'Color', col_cl61, 'LineWidth',1.5, ...
       'MarkerFaceColor', col_cl61, 'MarkerSize', 5, 'DisplayName','CL61');
loglog(tau_sec, allan_chm,  '-^', 'Color', col_chm,  'LineWidth',1.5, ...
       'MarkerFaceColor', col_chm,  'MarkerSize', 5, 'DisplayName','CHM15k');

i_anchor = find(~isnan(allan_cl61), 1);
if ~isempty(i_anchor)
    tau_ref = tau_sec(i_anchor:end);
    ref_line = allan_cl61(i_anchor) * sqrt(tau_sec(i_anchor) ./ tau_ref);
    loglog(tau_ref, ref_line, '--k', 'LineWidth', 1, 'DisplayName','\tau^{-1/2}');
end
hold off
grid on; box on
xlabel('Averaging duration \tau [s]')
ylabel('Std of \tau-averaged \beta_{att}, full range  [Mm^{-1} sr^{-1}]')
title('Spatial-spread metric')
legend('Location','best')
xticks(xt_sec); xticklabels(xt_label)
xlim([tau_sec(1) tau_sec(end)])

subplot(1,2,2)
loglog(tau_sec, adev_cl31, '-o', 'Color', col_cl31, 'LineWidth',1.5, ...
       'MarkerFaceColor', col_cl31, 'MarkerSize', 5, 'DisplayName','CL31');
hold on
loglog(tau_sec, adev_cl61, '-s', 'Color', col_cl61, 'LineWidth',1.5, ...
       'MarkerFaceColor', col_cl61, 'MarkerSize', 5, 'DisplayName','CL61');
loglog(tau_sec, adev_chm,  '-^', 'Color', col_chm,  'LineWidth',1.5, ...
       'MarkerFaceColor', col_chm,  'MarkerSize', 5, 'DisplayName','CHM15k');

i_anchor = find(~isnan(adev_cl61), 1);
if ~isempty(i_anchor)
    tau_ref = tau_sec(i_anchor:end);
    ref_line = adev_cl61(i_anchor) * sqrt(tau_sec(i_anchor) ./ tau_ref);
    loglog(tau_ref, ref_line, '--k', 'LineWidth', 1, 'DisplayName','\tau^{-1/2} (white)');
end

% Mark the drift knee per instrument.
local_mark_knee(tau_opt_31, adev_min_31, col_cl31);
local_mark_knee(tau_opt_61, adev_min_61, col_cl61);
local_mark_knee(tau_opt_ch, adev_min_ch, col_chm);
hold off
grid on; box on
xlabel('Averaging duration \tau [s]')
ylabel('Allan deviation of range-mean \beta_{att}  [Mm^{-1} sr^{-1}]')
title('Time-domain Allan deviation (knee = optimal averaging time)')
legend('Location','best')
xticks(xt_sec); xticklabels(xt_label)
xlim([tau_sec(1) tau_sec(end)])

sgtitle(['Noise vs averaging duration - ' tstr], 'FontWeight','bold')

%% Figure 8 Allan analysis resolved by altitude (6 tiles: overall + 5 levels)
% Error bars: chi-squared confidence interval on the OVERLAPPING Allan
% variance (IEEE Std 1139-2008; Riley, NIST SP 1065, 2008). The overlapping
% estimator maximises the equivalent degrees of freedom (edf); the white-FM
% edf (the tau^-1/2 branch) follows Howe, Allan & Barnes (1981). The variance
% is chi-squared with edf dof, so the 1-sigma interval on sigma(tau) is
%   [ sqrt(edf*avar/chi2inv(1-p/2,edf)), sqrt(edf*avar/chi2inv(p/2,edf)) ].
z_allan  = [500 1000 2000 3000 5000];   % altitudes [m]
half_win = 30;                          % +/- band [m]
conf_lvl = 0.683;                       % 1-sigma error bars (set 0.95 for 95%)
min_edf  = 4;                           % drop points with too few dof

panels(1) = struct('name','Overall (full range)', ...
    'ts31', mean(beta_cl31(:,mask_31),2,'omitnan'), ...
    'ts61', mean(beta_cl61(:,mask_61),2,'omitnan'), ...
    'tsch', mean(beta_chm( :,mask_ch),2,'omitnan'));
for j = 1:numel(z_allan)
    z0 = z_allan(j);
    panels(j+1) = struct('name', sprintf('%g m (\\pm%g m)', z0, half_win), ...
        'ts31', local_band_ts(beta_cl31, r_cl31, z0, half_win), ...
        'ts61', local_band_ts(beta_cl61, r_cl61, z0, half_win), ...
        'tsch', local_band_ts(beta_chm,  r_chm,  z0, half_win));
end

bigfig;
tl = tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');

for p = 1:numel(panels)
    nexttile

    [a31,l31,h31] = deal(nan(size(tau_sec)));
    [a61,l61,h61] = deal(nan(size(tau_sec)));
    [ach,lch,hch] = deal(nan(size(tau_sec)));
    for k = 1:numel(tau_sec)
        [a31(k),l31(k),h31(k)] = local_oadev_ci(panels(p).ts31, dt_cl31, tau_sec(k), conf_lvl, min_edf);
        [a61(k),l61(k),h61(k)] = local_oadev_ci(panels(p).ts61, dt_cl61, tau_sec(k), conf_lvl, min_edf);
        [ach(k),lch(k),hch(k)] = local_oadev_ci(panels(p).tsch, dt_chm,  tau_sec(k), conf_lvl, min_edf);
    end

    hold on
    errorbar(tau_sec, a31, a31-l31, h31-a31, '-o', 'Color',col_cl31, ...
        'MarkerFaceColor',col_cl31, 'MarkerSize',4, 'LineWidth',1.3, ...
        'CapSize',3, 'DisplayName','CL31');
    errorbar(tau_sec, a61, a61-l61, h61-a61, '-s', 'Color',col_cl61, ...
        'MarkerFaceColor',col_cl61, 'MarkerSize',4, 'LineWidth',1.3, ...
        'CapSize',3, 'DisplayName','CL61');
    errorbar(tau_sec, ach, ach-lch, hch-ach, '-^', 'Color',col_chm, ...
        'MarkerFaceColor',col_chm, 'MarkerSize',4, 'LineWidth',1.3, ...
        'CapSize',3, 'DisplayName','CHM15k');

    ia = find(isfinite(a61),1);
    if ~isempty(ia)
        tref = tau_sec(ia:end);
        plot(tref, a61(ia)*sqrt(tau_sec(ia)./tref), '--k', ...
            'LineWidth',0.8, 'DisplayName','\tau^{-1/2}');
    end
    hold off

    set(gca,'xscale','log','yscale','log')
    grid on; box on
    xticks(xt_sec); xticklabels(xt_label)
    xlim([tau_sec(1) tau_sec(end)])
    title(panels(p).name)
    if p == 1, legend('Location','best','FontSize',8); end
    if p > 3,  xlabel('Averaging duration \tau'); end
    if mod(p-1,3)==0, ylabel('Overlapping ADEV  [Mm^{-1} sr^{-1}]'); end
    ylim([1e-4 1e1])
end

title(tl, ['Altitude-resolved Allan deviation (\pm1\sigma \chi^2 CI) - ' tstr], ...
      'FontWeight','bold')
%% Optimal-averaging-time printout
fprintf('\n=== Time-domain Allan deviation: drift knee ===\n');
fprintf('%-8s  %-22s  %-16s\n', 'Instr', 'tau_opt (knee)', 'min Allan dev');
fprintf('%-8s  %-22s  %-16.4e\n', 'CL31', local_dur_str(tau_opt_31), adev_min_31);
fprintf('%-8s  %-22s  %-16.4e\n', 'CL61', local_dur_str(tau_opt_61), adev_min_61);
fprintf('%-8s  %-22s  %-16.4e\n', 'CHM15k', local_dur_str(tau_opt_ch), adev_min_ch);
fprintf('===============================================\n\n');

%% Non-range-corrected signal
r2_cl31 = r_cl31'.^2;
r2_cl61 = r_cl61'.^2;
r2_chm  = r_chm'.^2;

praw_cl31 = beta_cl31 ./ r2_cl31;
praw_cl61 = beta_cl61 ./ r2_cl61;
praw_chm  = beta_chm  ./ r2_chm;

med_praw_cl31 = median(praw_cl31, 1, 'omitnan');
std_praw_cl31 = std(   praw_cl31, 0, 1, 'omitnan');
rms_praw_cl31 = sqrt(mean(praw_cl31.^2, 1, 'omitnan'));

med_praw_cl61 = median(praw_cl61, 1, 'omitnan');
std_praw_cl61 = std(   praw_cl61, 0, 1, 'omitnan');
rms_praw_cl61 = sqrt(mean(praw_cl61.^2, 1, 'omitnan'));

med_praw_chm  = median(praw_chm,  1, 'omitnan');
std_praw_chm  = std(   praw_chm,  0, 1, 'omitnan');
rms_praw_chm  = sqrt(mean(praw_chm.^2, 1, 'omitnan'));

%% Figure 8 : raw signal profiles
bigfig;

subplot(1,3,1)
hold on
plot(med_praw_cl31, r_cl31, '-', 'Color', col_cl31, 'LineWidth', 1.5, 'DisplayName','CL31');
plot(med_praw_cl61, r_cl61, '-', 'Color', col_cl61, 'LineWidth', 1.5, 'DisplayName','CL61');
plot(med_praw_chm,  r_chm,  '-', 'Color', col_chm,  'LineWidth', 1.5, 'DisplayName','CHM15k');
xline(0, 'k:', 'HandleVisibility','off')
hold off
box on; grid on
xlabel('Median raw signal  [Mm^{-1} sr^{-1} m^{-2}]')
ylabel('Range [m]')
title('Median of non-range-corrected signal')
ylim([0 range_top])
legend('Location','best')

subplot(1,3,2)
hold on
plot(abs(std_praw_cl31), r_cl31, '-', 'Color', col_cl31, 'LineWidth', 1.5, 'DisplayName','CL31');
plot(abs(std_praw_cl61), r_cl61, '-', 'Color', col_cl61, 'LineWidth', 1.5, 'DisplayName','CL61');
plot(abs(std_praw_chm),  r_chm,  '-', 'Color', col_chm,  'LineWidth', 1.5, 'DisplayName','CHM15k');
hold off
box on; grid on
set(gca,'xscale','log')
xlabel('\sigma raw signal  [Mm^{-1} sr^{-1} m^{-2}]')
ylabel('Range [m]')
title('Std of non-range-corrected signal')
ylim([0 range_top])
legend('Location','best')

subplot(1,3,3)
hold on
plot(rms_praw_cl31, r_cl31, '-', 'Color', col_cl31, 'LineWidth', 1.5, 'DisplayName','CL31');
plot(rms_praw_cl61, r_cl61, '-', 'Color', col_cl61, 'LineWidth', 1.5, 'DisplayName','CL61');
plot(rms_praw_chm,  r_chm,  '-', 'Color', col_chm,  'LineWidth', 1.5, 'DisplayName','CHM15k');
hold off
box on; grid on
set(gca,'xscale','log')
xlabel('RMS raw signal  [Mm^{-1} sr^{-1} m^{-2}]')
ylabel('Range [m]')
title('RMS of non-range-corrected signal')
ylim([0 range_top])
legend('Location','best')

sgtitle(['Non-range-corrected signal - ' tstr], 'FontWeight','bold')

%% DC bias convergence
t_meas_max = seconds(max([cl31.dtime(end)-cl31.dtime(1), ...
                          cl61.dtime(end)-cl61.dtime(1), ...
                          chm.dtime( end)-chm.dtime( 1)]));

t_meas = unique(round(logspace(log10(60), log10(t_meas_max), 30)));

run_med_cl31 = nan(size(t_meas));
run_med_cl61 = nan(size(t_meas));
run_med_chm  = nan(size(t_meas));

t0_cl31 = cl31.dtime(1);
t0_cl61 = cl61.dtime(1);
t0_chm  = chm.dtime(1);

for k = 1:numel(t_meas)
    tcut = seconds(t_meas(k));

    sel = cl31.dtime - t0_cl31 <= tcut;
    if any(sel)
        block = praw_cl31(sel, :);
        run_med_cl31(k) = median(block(:), 'omitnan');
    end

    sel = cl61.dtime - t0_cl61 <= tcut;
    if any(sel)
        block = praw_cl61(sel, :);
        run_med_cl61(k) = median(block(:), 'omitnan');
    end

    sel = chm.dtime - t0_chm <= tcut;
    if any(sel)
        block = praw_chm(sel, :);
        run_med_chm(k) = median(block(:), 'omitnan');
    end
end

%% Figure 9 : DC bias convergence
bigfig;

subplot(1,2,1)
semilogx(t_meas, run_med_cl31, '-o', 'Color', col_cl31, 'LineWidth',1.5, ...
         'MarkerFaceColor', col_cl31, 'MarkerSize', 5, 'DisplayName','CL31');
hold on
semilogx(t_meas, run_med_cl61, '-s', 'Color', col_cl61, 'LineWidth',1.5, ...
         'MarkerFaceColor', col_cl61, 'MarkerSize', 5, 'DisplayName','CL61');
semilogx(t_meas, run_med_chm,  '-^', 'Color', col_chm,  'LineWidth',1.5, ...
         'MarkerFaceColor', col_chm,  'MarkerSize', 5, 'DisplayName','CHM15k');
yline(0, 'k:', 'HandleVisibility','off')
hold off
grid on; box on
xlabel('Measurement duration T [s]')
ylabel('Running median of raw signal  [Mm^{-1} sr^{-1} m^{-2}]')
title('DC bias estimate vs measurement time')
xticks(xt_sec); xticklabels(xt_label)
xlim([t_meas(1) t_meas(end)])
legend('Location','best')

subplot(1,2,2)
final_31 = run_med_cl31(end); final_61 = run_med_cl61(end); final_ch = run_med_chm(end);
err_31 = abs(run_med_cl31 - final_31);
err_61 = abs(run_med_cl61 - final_61);
err_ch = abs(run_med_chm  - final_ch );
loglog(t_meas, err_31, '-o', 'Color', col_cl31, 'LineWidth',1.5, ...
       'MarkerFaceColor', col_cl31, 'MarkerSize', 5, 'DisplayName','CL31');
hold on
loglog(t_meas, err_61, '-s', 'Color', col_cl61, 'LineWidth',1.5, ...
       'MarkerFaceColor', col_cl61, 'MarkerSize', 5, 'DisplayName','CL61');
loglog(t_meas, err_ch, '-^', 'Color', col_chm,  'LineWidth',1.5, ...
       'MarkerFaceColor', col_chm,  'MarkerSize', 5, 'DisplayName','CHM15k');
hold off
grid on; box on
xlabel('Measurement duration T [s]')
ylabel('|median(T) - median(T_{max})|  [Mm^{-1} sr^{-1} m^{-2}]')
title('Deviation from final value (proxy for uncertainty)')
xticks(xt_sec); xticklabels(xt_label)
xlim([t_meas(1) t_meas(end)])
legend('Location','best')

sgtitle(['DC bias of dark signal vs measurement duration - ' tstr], 'FontWeight','bold')

%% Required-time printout
target_abs = 1e-7;

fprintf('\n=== Measurement time to determine DC bias to %.1e ===\n', target_abs);
fprintf('%-8s  %-22s  %-14s\n', 'Instr', 'T_required', 'final median');
fprintf('%-8s  %-22s  %-14.3e\n', 'CL31', ...
        local_time_to_stable(t_meas, err_31, target_abs), final_31);
fprintf('%-8s  %-22s  %-14.3e\n', 'CL61', ...
        local_time_to_stable(t_meas, err_61, target_abs), final_61);
fprintf('%-8s  %-22s  %-14.3e\n', 'CHM15k', ...
        local_time_to_stable(t_meas, err_ch, target_abs), final_ch);
fprintf('=========================================================\n\n');

%% Aerosol detection threshold from measured dark noise
% Hypothesis: the dark-measurement noise is the limiting factor. The
% per-gate temporal std of attenuated backscatter measured with covered
% telescopes (std_cl31/std_cl61/std_chm, in Mm^-1 sr^-1) is taken as the
% noise floor sigma_dark(r). The minimum DETECTABLE attenuated backscatter
% at range r for a given SNR is simply:
%       beta_att_min(r) = SNR * sigma_dark(r).
%
% Averaging: sigma scales as tau^-1/2 (white-noise regime, confirmed by the
% Allan analysis below the drift knee at ~2 h). We reference sigma_dark to
% its native sampling time and scale to any requested averaging time tau.
%
% From attenuated to true aerosol backscatter: beta_att = beta * T_mol^2,
% where T_mol^2 = exp(-2 * integral_0^r alpha_mol dr') is the two-way
% molecular transmission. To recover the true aerosol backscatter the
% instrument can detect, we divide the attenuated threshold by T_mol^2
% (aerosol two-way transmission neglected at threshold: thin layer limit).
%
% Extinction:  alpha_min(r) = LR * beta_min(r)        [LR = lidar ratio, sr]
% Mass conc.:  M_min(r)     = alpha_min(r) / MEC      [MEC = mass extinction
%                                                      coefficient, m^2/g]

%% Detection-threshold parameters (edit here)
SNR_main   = 3;             % SNR detection threshold (3-sigma)
SNR_list   = [2 3 5];       % SNR levels for the sensitivity panel
TAU_LIST   = [300 1800];    % averaging times [s] shown side by side (5min, 30min)

% Volcanic aerosol optical/microphysical assumptions.
% Scenario A: fine volcanic ash, classic EARLINET values
%   (Ansmann et al. 2011; Eyjafjallajokull 2010). LR at 532/910-1064 nm.
% Scenario B: alternative (coarser / different composition).
volc(1).name = 'Ash: LR=60 sr, MEC=0.60 m^2/g';
volc(1).LR   = 60;          % sr
volc(1).MEC  = 0.60;        % m^2/g
volc(2).name = 'Ash: LR=50 sr, MEC=0.80 m^2/g';
volc(2).LR   = 50;
volc(2).MEC  = 0.80;
volc_main = 1;              % scenario used for the headline mass figure

% ICAO volcanic-ash contamination zones (mass concentration thresholds).
% Source: ICAO EUR/NAT Volcanic Ash Contingency Plan, ICAO EUR Doc 019 /
% NAT Doc 006 Part II (Dec 2010 edition), as used by the London VAAC ash
% concentration charts:
%   low contamination     :  0.2e-3 .. 2e-3  g/m^3  (200 .. 2000 ug/m^3)
%   medium contamination  :  2e-3   .. 4e-3  g/m^3  (2000 .. 4000 ug/m^3)
%   high contamination     :  > 4e-3          g/m^3  (> 4000 ug/m^3)
% So: 200 = outer edge of the low zone, 2000 = upper bound of low /
% start of medium, 4000 = lower bound of high. (The newer QVA framework,
% ICAO Annex 3 Amdt 82 applicable Nov 2025, uses 0.2/2/5/10 mg/m^3 instead.)
icao.level = [200 2000 4000];   % ug/m^3
icao.label = {'200 \mug/m^3 (low edge)', ...
              '2000 \mug/m^3 (low/medium)', ...
              '4000 \mug/m^3 (medium/high)'};
icao.style = {'--', '-.', ':'};

% Standard-atmosphere molecular extinction at the ceilometer wavelengths.
% alpha_mol(z) ~ alpha_mol(0) * exp(-z/H_scale), Rayleigh.
% alpha_mol(0) from molecular backscatter * (8*pi/3); values are typical
% mid-latitude clear-air figures and only weakly affect the threshold.
H_scale       = 8000;       % atmospheric scale height [m]
alpha_mol0_910  = 1.0e-6;   % m^-1 at 910 nm  (CL31/CL61)
alpha_mol0_1064 = 5.0e-7;   % m^-1 at 1064 nm (CHM15k)

%% Native per-gate dark noise (scaled to each tau at plot time)
% std_* were computed over the full record at native resolution, so they
% represent the single-shot (native dt) noise per gate. The tau^-1/2
% white-noise scaling (valid below the ~2 h Allan drift knee) is applied
% per averaging time when plotting:  sigma(tau) = sigma0 * sqrt(dt/tau).
sig_dark_31 = std_cl31(:);   % Mm^-1 sr^-1, per gate, native dt
sig_dark_61 = std_cl61(:);
sig_dark_ch = std_chm(:);

%% Two-way molecular transmission per instrument
Tmol2_31 = exp(-2 * alpha_mol0_910  * H_scale * (1 - exp(-r_cl31/H_scale)));
Tmol2_61 = exp(-2 * alpha_mol0_910  * H_scale * (1 - exp(-r_cl61/H_scale)));
Tmol2_ch = exp(-2 * alpha_mol0_1064 * H_scale * (1 - exp(-r_chm /H_scale)));

%% Figure 10 : minimum detectable backscatter vs altitude
bigfig;

inst(1) = struct('r',r_cl31,'sig',sig_dark_31,'dt',dt_cl31,'Tmol2',Tmol2_31, ...
                 'col',col_cl31,'name','CL31');
inst(2) = struct('r',r_cl61,'sig',sig_dark_61,'dt',dt_cl61,'Tmol2',Tmol2_61, ...
                 'col',col_cl61,'name','CL61');
inst(3) = struct('r',r_chm, 'sig',sig_dark_ch,'dt',dt_chm, 'Tmol2',Tmol2_ch, ...
                 'col',col_chm, 'name','CHM15k');

% One panel per averaging time, all three instruments at SNR_main.
% A thin grey line shows the SNR=2 / SNR=5 envelope for CL61 as a guide.
ntau = numel(TAU_LIST);
for p = 1:ntau
    tau = TAU_LIST(p);
    subplot(1, ntau, p)
    hold on
    for i = 1:numel(inst)
        sig_tau  = inst(i).sig .* sqrt(inst(i).dt ./ max(tau, inst(i).dt));
        beta_min = SNR_main * sig_tau ./ inst(i).Tmol2;   % true aerosol beta
        plot(beta_min, inst(i).r, '-', 'Color', inst(i).col, 'LineWidth', 2.0, ...
             'DisplayName', inst(i).name);
    end
    hold off
    box on; grid on
    set(gca,'xscale','log')
    xlabel('\beta_{min} aerosol  [Mm^{-1} sr^{-1}]')
    ylabel('Altitude [m]')
    title(sprintf('SNR=%d, \\tau = %s', SNR_main, local_human_dur(tau)))
    ylim([0 range_top])
    xlim([1e-3 1e1])
    legend('Location','best')
end

sgtitle(['Aerosol detection threshold - backscatter - ' tstr], 'FontWeight','bold')

%% Figure 11 : minimum detectable extinction and mass (volcanic)
bigfig;

LR  = volc(volc_main).LR;
MEC = volc(volc_main).MEC;

% Line style per averaging time: solid = longest tau, dashed = shorter.
[~, ord] = sort(TAU_LIST, 'descend');
tau_styles = {'-', '--', ':'};

% Extinction threshold: alpha_min = LR * beta_min.
% beta_min in Mm^-1 sr^-1; LR in sr gives alpha in Mm^-1.
subplot(1,3,1)
hold on
for i = 1:numel(inst)
    for q = 1:numel(ord)
        tau = TAU_LIST(ord(q));
        sig_tau   = inst(i).sig .* sqrt(inst(i).dt ./ max(tau, inst(i).dt));
        beta_min  = SNR_main * sig_tau ./ inst(i).Tmol2;     % Mm^-1 sr^-1
        alpha_min = LR * beta_min;                            % Mm^-1
        if q == 1, hv = 'on'; else, hv = 'off'; end
        plot(alpha_min, inst(i).r, tau_styles{q}, 'Color', inst(i).col, ...
             'LineWidth', 2.0, 'HandleVisibility', hv, 'DisplayName', inst(i).name);
    end
end
hold off
box on; grid on
set(gca,'xscale','log')
xlabel('\alpha_{min}  [Mm^{-1}]')
ylabel('Altitude [m]')
title(sprintf('Min detectable extinction (LR=%g sr)', LR))
ylim([0 range_top])
xlim([1e-2 1e3])
legend('Location','best')

% Mass concentration threshold: M = alpha / MEC.
% alpha in Mm^-1 = 1e-6 m^-1; MEC in m^2/g; M in g/m^3 -> *1e6 -> ug/m^3.
subplot(1,3,2)
hold on
for i = 1:numel(inst)
    for q = 1:numel(ord)
        tau = TAU_LIST(ord(q));
        sig_tau   = inst(i).sig .* sqrt(inst(i).dt ./ max(tau, inst(i).dt));
        beta_min  = SNR_main * sig_tau ./ inst(i).Tmol2;     % Mm^-1 sr^-1
        alpha_min = LR * beta_min * 1e-6;                     % m^-1
        M_min_ug  = (alpha_min / MEC) * 1e6;                  % ug/m^3
        if q == 1, hv = 'on'; else, hv = 'off'; end
        plot(M_min_ug, inst(i).r, tau_styles{q}, 'Color', inst(i).col, ...
             'LineWidth', 2.0, 'HandleVisibility', hv, 'DisplayName', inst(i).name);
    end
end
% ICAO contamination-zone thresholds (see definition block above).
for g = 1:numel(icao.level)
    xline(icao.level(g), ['k' icao.style{g}], 'LineWidth', 1, ...
          'DisplayName', icao.label{g});
end
hold off
box on; grid on
set(gca,'xscale','log')
xlabel('M_{min}  [\mug/m^3]')
ylabel('Altitude [m]')
title(sprintf('Min detectable mass (MEC=%g m^2/g; solid=%s, dash=%s)', MEC, ...
      local_human_dur(TAU_LIST(ord(1))), local_human_dur(TAU_LIST(ord(end)))))
ylim([0 range_top])
legend('Location','best')
xlim([1e-2 1e3])

% Mass threshold: volcanic scenario sensitivity, CL61 only, longest tau.
subplot(1,3,3)
tau_sens = TAU_LIST(ord(1));
sig_tau_61 = sig_dark_61 .* sqrt(dt_cl61 ./ max(tau_sens, dt_cl61));
hold on
for v = 1:numel(volc)
    beta_min  = SNR_main * sig_tau_61 ./ Tmol2_61;           % Mm^-1 sr^-1
    alpha_min = volc(v).LR * beta_min * 1e-6;                 % m^-1
    M_min_ug  = (alpha_min / volc(v).MEC) * 1e6;              % ug/m^3
    plot(M_min_ug, r_cl61, '-', 'LineWidth', 1.8, 'DisplayName', volc(v).name);
end
for g = 1:numel(icao.level)
    xline(icao.level(g), ['k' icao.style{g}], 'LineWidth', 1, ...
          'HandleVisibility', 'off');
end
hold off
box on; grid on
set(gca,'xscale','log')
xlabel('M_{min}  [\mug/m^3]')
ylabel('Altitude [m]')
title(sprintf('CL61: volcanic scenario sensitivity (\\tau=%s)', ...
      local_human_dur(tau_sens)))
ylim([0 range_top])
legend('Location','best')

sgtitle({['Aerosol detection threshold - extinction & mass (SNR=' ...
          num2str(SNR_main) ') - ' tstr], ...
         'ICAO ash zones: EUR Doc 019 / NAT Doc 006 (2010) - 200/2000/4000 \mug/m^3'}, ...
         'FontWeight','bold')

%% Detection-threshold printout (at selected altitudes)
z_probe = [500 1000 2000 3000 5000];
for tau = TAU_LIST
    fprintf('\n=== Aerosol detection thresholds, SNR=%d, tau=%s, %s ===\n', ...
            SNR_main, local_human_dur(tau), volc(volc_main).name);
    fprintf('(beta in Mm^-1 sr^-1, alpha in Mm^-1, M in ug/m^3)\n');
    fprintf('%-8s  %-7s  %-12s  %-12s  %-12s\n', ...
            'Instr','z[m]','beta_min','alpha_min','M_min');
    for i = 1:numel(inst)
        for z = z_probe
            [~, iz] = min(abs(inst(i).r - z));
            sig_tau = inst(i).sig(iz) * sqrt(inst(i).dt / max(tau, inst(i).dt));
            bmin = SNR_main * sig_tau / inst(i).Tmol2(iz);
            amin = LR * bmin;
            mmin = (LR * bmin * 1e-6 / MEC) * 1e6;
            fprintf('%-8s  %-7d  %-12.4e  %-12.4e  %-12.2f\n', ...
                    inst(i).name, z, bmin, amin, mmin);
        end
    end
    fprintf('=====================================================\n');
end
fprintf('\n');

%% Save all figures
out_dir = fullfile(pwd, 'figs_ceilo_dark');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
tag = sprintf('%s_%s', datestr(t1,'yyyymmdd_HHMM'), datestr(t2,'HHMM'));

figs = findall(0, 'Type','figure');
figs = flipud(figs);
names = {'01_timeheight', '02_profile_noise_rms', '03_nearrange', ...
         '04_noise_vs_T',  '05_noise_T_timeseries', '06_cl61_noise_check', ...
         '07_allan_beta',  '08_raw_profile',        '09_raw_median_tau', ...
         '10_detect_backscatter', '11_detect_extinction_mass'};

for k = 1:numel(figs)
    fh = figs(k);
    if k <= numel(names)
        fname = sprintf('%s_%s', names{k}, tag);
    else
        fname = sprintf('fig%02d_%s', k, tag);
    end
    set(fh, 'Color', 'w');
    try
        exportgraphics(fh, fullfile(out_dir, [fname '.png']), 'Resolution', 150);
    catch
        saveas(fh, fullfile(out_dir, [fname '.png']));
    end
end
fprintf('Saved %d figures to %s\n', numel(figs), out_dir);



%% local helpers
function r = local_corr(x, y)
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    if nnz(ok) < 10
        r = NaN; return;
    end
    R = corrcoef(x(ok), y(ok));
    r = R(1,2);
end

function [allan_std, nwin] = local_allan_noise(beta_block, dt, tau)
    nwin_samples = max(1, round(tau / dt));
    Nt = size(beta_block, 1);
    Nseg = floor(Nt / nwin_samples);

    if Nseg < 2
        allan_std = NaN; nwin = Nseg; return;
    end

    trimmed = beta_block(1:Nseg*nwin_samples, :);
    Nrange = size(trimmed, 2);
    rs = reshape(trimmed, [nwin_samples, Nseg, Nrange]);
    win_mean = squeeze(mean(rs, 1, 'omitnan'));
    s_per_seg = std(win_mean, 0, 2, 'omitnan');
    allan_std = mean(s_per_seg, 'omitnan');
    nwin = Nseg;
end

function [adev, nseg] = local_allan_time(x, dt, tau)
    % Non-overlapping Allan deviation of a 1-D time series x sampled at dt.
    % adev(tau) = sqrt( 0.5 * mean( diff(segment_means)^2 ) ).
    x = x(:);
    m = max(1, round(tau / dt));
    nseg = floor(numel(x) / m);
    if nseg < 2
        adev = NaN; return;
    end
    x = x(1:nseg*m);
    seg_means = mean(reshape(x, [m, nseg]), 1, 'omitnan');
    d = diff(seg_means);
    adev = sqrt(0.5 * mean(d.^2, 'omitnan'));
end

function [tau_opt, adev_min] = local_allan_min(tau_sec, adev)
    ok = isfinite(adev);
    if ~any(ok)
        tau_opt = NaN; adev_min = NaN; return;
    end
    a = adev; a(~ok) = Inf;
    [adev_min, idx] = min(a);
    tau_opt = tau_sec(idx);
end

function local_mark_knee(tau_opt, adev_min, col)
    if isnan(tau_opt) || isnan(adev_min)
        return;
    end
    loglog(tau_opt, adev_min, 'p', 'MarkerSize', 14, ...
           'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
           'LineWidth', 1, 'HandleVisibility', 'off');
end

function s = local_dur_str(sec)
    if isnan(sec)
        s = 'n/a'; return;
    end
    s = sprintf('%.0fs (%s)', sec, local_human_dur(sec));
end

function h = local_human_dur(s)
    if s < 60
        h = sprintf('%.0fs', s);
    elseif s < 3600
        h = sprintf('%.1fmin', s/60);
    else
        h = sprintf('%.1fh', s/3600);
    end
end

function s = local_time_to_stable(t_meas, err, target)
    ok = ~isnan(err);
    if ~any(ok)
        s = 'n/a'; return;
    end
    below = err <= target;
    last_bad = find(~below & ok, 1, 'last');
    if isempty(last_bad)
        first_ok = find(ok, 1, 'first');
        s = sprintf('%.0fs (%s)', t_meas(first_ok), local_human_dur(t_meas(first_ok)));
        return;
    end
    if last_bad >= numel(t_meas)
        s = sprintf('>%.0fs (need longer)', t_meas(end));
        return;
    end
    s = sprintf('%.0fs (%s)', t_meas(last_bad+1), ...
                local_human_dur(t_meas(last_bad+1)));
end

function ts = local_band_ts(beta, r, z0, halfwin)
    sel = r >= (z0 - halfwin) & r <= (z0 + halfwin);
    if ~any(sel)
        ts = nan(size(beta,1),1);
    else
        ts = mean(beta(:, sel), 2, 'omitnan');
    end
end

function [adev, lo, hi] = local_oadev_ci(y, dt, tau, conf, min_edf)
    % Overlapping Allan deviation with chi-squared confidence interval.
    % Refs: Allan & Barnes (1981); IEEE Std 1139-2008; Riley, NIST SP 1065
    % (2008). White-FM edf from Howe, Allan & Barnes (1981).
    adev = NaN; lo = NaN; hi = NaN;
    y = fillmissing(y(:), 'linear', 'EndValues','none');
    ok = isfinite(y);
    if mean(ok) < 0.5, return; end
    y = y(ok);
    M = numel(y);                      % frequency-type samples
    m = max(1, round(tau/dt));         % averaging factor
    if M < 2*m + 1, return; end
    x   = [0; cumsum(y)] * dt;         % integrated (phase), M+1 points
    idx = 1:(numel(x) - 2*m);
    d   = x(idx + 2*m) - 2*x(idx + m) + x(idx);
    avar = sum(d.^2) / (2 * tau^2 * numel(idx));
    adev = sqrt(avar);
    edf = (3*(M-1)/(2*m) - 2*(M-2)/M) * (4*m^2) / (4*m^2 + 5);   % white FM
    if ~isfinite(edf) || edf < min_edf, return; end
    alpha = 1 - conf;
    hi = sqrt(edf * avar / chi2inv(alpha/2,     edf));   % upper sigma bound
    lo = sqrt(edf * avar / chi2inv(1 - alpha/2, edf));   % lower sigma bound
end