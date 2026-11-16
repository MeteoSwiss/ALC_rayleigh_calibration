%% export_11yr_archive.m
% Export the 11-year daily CHM15k housekeeping + lidar-constant tables (A:\CHM15k\<stn>\,
% 2013-2024, the dataset behind analyze_chm15k_stability.m) to CSV so the Python report-19
% analysis can control the temperature dependence for season. One CSV per station.

outdir = 'C:\DATA\Projects\202610_CHM15k_stability\archive_11yr';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
stations = {'pay','kse','eri','ber','dav','cbw','gal','ino','joy','lei','lin','mhd','sir'};
for s = 1:numel(stations)
    f = fullfile('A:\CHM15k', stations{s}, ['chm_' stations{s} '_1d_01-Jan-2013-31-Dec-2024.mat']);
    if ~exist(f, 'file')
        fprintf('%s: no file\n', stations{s});
        continue;
    end
    S = load(f, 'T');
    T = S.T;
    if istimetable(T)
        T = timetable2table(T);
    end
    % keep only numeric / datetime columns so the CSV stays flat
    keep = false(1, width(T));
    for k = 1:width(T)
        v = T.(T.Properties.VariableNames{k});
        keep(k) = (isnumeric(v) || islogical(v) || isdatetime(v)) && size(v, 2) == 1;
    end
    T = T(:, keep);
    writetable(T, fullfile(outdir, ['chm_' stations{s} '_1d.csv']));
    fprintf('%s: %d rows, %d cols\n', stations{s}, height(T), width(T));
end
fprintf('EXPORT_DONE\n');
