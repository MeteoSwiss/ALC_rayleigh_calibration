% Export per-instrument MATLAB Rayleigh results (daily_dates, daily_C, daily_C_std)
% to CSV so the Python side can do an exact day-by-day comparison. The per-station
% rayleigh_<WMO>_<id>.mat stores daily_dates as a MATLAB datetime (an MCOS object
% scipy cannot decode), so the conversion to yyyymmdd strings is done here in MATLAB.

src = 'A:\E-PROFILE_L2_Calibration\rayleigh_per_station';
out = 'D:\E-PROFILE_calibration_rayleigh\matlab_daily_export';
if ~isfolder(out)
    mkdir(out);
end

files = dir(fullfile(src, 'rayleigh_*.mat'));
nok = 0;
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    S = load(f);
    if ~isfield(S, 'daily_dates') || ~isfield(S, 'daily_C')
        continue;
    end
    d = S.daily_dates(:);
    C = S.daily_C(:);
    if isfield(S, 'daily_C_std')
        Cstd = S.daily_C_std(:);
    else
        Cstd = nan(size(C));
    end

    % Keep only valid (non-NaT) dates; align lengths defensively.
    n = min([numel(d), numel(C), numel(Cstd)]);
    d = d(1:n);
    C = C(1:n);
    Cstd = Cstd(1:n);
    valid = ~isnat(d);
    d = d(valid);
    C = C(valid);
    Cstd = Cstd(valid);

    dateStr = string(datestr(d, 'yyyymmdd'));
    T = table(dateStr, C, Cstd, 'VariableNames', {'date', 'C', 'C_std'});

    [~, stem, ~] = fileparts(files(k).name);   % rayleigh_<WMO>_<id>
    writetable(T, fullfile(out, [stem '.csv']));
    nok = nok + 1;
end
fprintf('exported %d / %d station files to %s\n', nok, numel(files), out);
