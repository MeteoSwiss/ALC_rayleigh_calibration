% ===== Provenance (paper code archive) =====
% Paper Figure 1 (right): number of E-Profile stations with valid L2 data per month, by instrument type.
% Run in MATLAB with no arguments from the output folder; writes station_counts_by_month.png.
% Input: A:\E-PROFILE_L2_monthly\<wigos>\**\L2_*_<id><YYYYMM>.nc (type from global attribute instrument_type),
% Nov 2016 - Dec 2025. The paper image is a cropped screenshot of this figure (made 2026-02-16).
% Source: MATLAB/ALC repository (service.meteoswiss.ch/git/remotesensing/alc.git), liquid_cloud_calibration/,
% commit 89de53d. Copied here unchanged on 2026-10-20 so that the paper code lives in one repository.
% ==========================================
%% count_stations_by_month.m
% Script to count the number of stations with valid data every month
% and display cumulative bar charts colored by instrument type
%
% This script:
% 1. Scans all station folders in the data directory
% 2. Detects instrument type from NetCDF file attributes
% 3. Counts stations with valid data for each month
% 4. Creates cumulative stacked bar charts by instrument type

%% Clear workspace
clear; close all; clc;

%% CONFIGURATION

% Base data directory containing station folders
base_data_dir = 'A:\E-PROFILE_L2_monthly\';

% Date range to process
start_date = datetime(2016, 11, 1);
end_date = datetime(2025, 12, 31);

% Colors for known instrument types (others will get auto-assigned colors)
instrument_colors = struct();
instrument_colors.CL31 = [0.2 0.4 0.8];      % Blue
instrument_colors.CL51 = [0.8 0.4 0.2];      % Orange  
instrument_colors.CL61 = [0.2 0.7 0.3];      % Green
instrument_colors.CHM15k = [0.6 0.2 0.6];    % Purple
instrument_colors.CHM15KX = [0.8 0.2 0.8];   % Magenta
instrument_colors.MiniMPL = [0.9 0.6 0.1];   % Gold
instrument_colors.CS135 = [0.1 0.6 0.6];     % Teal

%% MAIN EXECUTION

fprintf('=== Station Count by Month ===\n\n');
fprintf('Scanning for stations in: %s\n', base_data_dir);

% Generate month range
d0 = dateshift(start_date, 'start', 'month');
d1 = dateshift(end_date, 'start', 'month');
month_dates = (d0:calmonths(1):d1).';
n_months = length(month_dates);

fprintf('Date range: %s to %s (%d months)\n\n', ...
    datestr(start_date, 'yyyy-mm'), datestr(end_date, 'yyyy-mm'), n_months);

% List all station folders
station_dirs = dir(base_data_dir);
station_dirs = station_dirs([station_dirs.isdir]);
station_dirs = station_dirs(~ismember({station_dirs.name}, {'.', '..'}));

% Filter for valid WIGOS station IDs (folders starting with 0-)
valid_stations = {};
for i = 1:length(station_dirs)
    name = station_dirs(i).name;
    if startsWith(name, '0-')
        valid_stations{end+1} = name;
    end
end

fprintf('Found %d station folders\n', length(valid_stations));

% Scan all instruments across all stations (discovers all instrument types)
fprintf('Scanning for instruments...\n');
all_instruments = scan_all_instruments(base_data_dir, valid_stations);

% Get unique instrument types found
instrument_types = unique({all_instruments.instrument_type});
instrument_types = instrument_types(~cellfun(@isempty, instrument_types));

fprintf('Found %d instruments across all stations\n', length(all_instruments));
fprintf('Instrument types detected: %s\n\n', strjoin(instrument_types, ', '));

% Count stations per month per instrument type
[monthly_counts, instrument_types] = count_stations_per_month(...
    all_instruments, month_dates, instrument_types);

%% Display summary table
display_summary(month_dates, monthly_counts, instrument_types);

%% Create cumulative bar chart
plot_cumulative_bars(month_dates, monthly_counts, instrument_types, instrument_colors);

%% FUNCTIONS

function all_instruments = scan_all_instruments(base_data_dir, valid_stations)
% SCAN_ALL_INSTRUMENTS - Find all instruments at all stations
% Returns a struct array with instrument info and available months
% Automatically discovers all instrument types

all_instruments = struct('wmo', {}, 'instrument_id', {}, 'instrument_type', {}, ...
                         'available_months', {}, 'lat', {}, 'lon', {}, 'station_name', {});

total_instruments_found = 0;
total_instruments_recognized = 0;

for s = 1:length(valid_stations)
    wmo_code = valid_stations{s};
    station_path = fullfile(base_data_dir, wmo_code);
    
    [instruments, n_found, n_recognized] = detect_instruments_with_months(station_path, wmo_code);
    total_instruments_found = total_instruments_found + n_found;
    total_instruments_recognized = total_instruments_recognized + n_recognized;
    
    for i = 1:length(instruments)
        all_instruments(end+1) = instruments(i);
    end
    
    % Progress indicator
    if mod(s, 50) == 0
        fprintf('  Scanned %d/%d stations...\n', s, length(valid_stations));
    end
end

fprintf('  Scanned %d/%d stations (complete)\n', length(valid_stations), length(valid_stations));
fprintf('  Total instrument IDs found in files: %d\n', total_instruments_found);
fprintf('  Successfully recognized instruments: %d\n', total_instruments_recognized);
if total_instruments_found > total_instruments_recognized
    fprintf('  WARNING: %d instruments could not be classified (see warnings above)\n', ...
            total_instruments_found - total_instruments_recognized);
end

end

function [instruments, n_found, n_recognized] = detect_instruments_with_months(station_path, wmo_code)
% DETECT_INSTRUMENTS_WITH_MONTHS - Find all instruments at a station with their available months
% Automatically detects instrument type from NetCDF attributes
% Returns:
%   instruments - struct array of recognized instruments
%   n_found - total number of instrument IDs found in filenames
%   n_recognized - number of instruments successfully classified

instruments = struct('wmo', {}, 'instrument_id', {}, 'instrument_type', {}, ...
                     'available_months', {}, 'lat', {}, 'lon', {}, 'station_name', {});
n_found = 0;
n_recognized = 0;

% Find all nc files across all year folders
all_files = dir(fullfile(station_path, '**', '*.nc'));

if isempty(all_files)
    return;
end

% Build list of unique instrument IDs from filenames with their months
% Filename pattern: L2_0-xxx-xxx-xxx_Y202505.nc
instrument_files = struct('instrument_id', {}, 'files', {}, 'months', {});
found_ids = {};

for f = 1:length(all_files)
    [~, fname, ~] = fileparts(all_files(f).name);
    
    % Parse filename: L2_0-xxx-xxx-xxx_Y202505
    parts = strsplit(fname, '_');
    if length(parts) < 3
        continue;
    end
    
    % Last part should be instrument_id + YYYYMM (e.g., A202505)
    last_part = parts{end};
    if length(last_part) < 7
        continue;
    end
    
    % First character is instrument ID, rest is YYYYMM
    instrument_id = last_part(1);
    yyyymm = last_part(2:end);
    
    % Validate: instrument_id should be alphanumeric, yyyymm should be 6 digits
    if ~isstrprop(instrument_id, 'alphanum') || length(yyyymm) ~= 6 || ~all(isstrprop(yyyymm, 'digit'))
        continue;
    end
    
    % Parse month
    yr = str2double(yyyymm(1:4));
    mo = str2double(yyyymm(5:6));
    month_dt = datetime(yr, mo, 1);
    
    % Track this file for this instrument
    idx = find(strcmp(found_ids, instrument_id), 1);
    if isempty(idx)
        found_ids{end+1} = instrument_id;
        instrument_files(end+1).instrument_id = instrument_id;
        instrument_files(end).files = {fullfile(all_files(f).folder, all_files(f).name)};
        instrument_files(end).months = month_dt;
    else
        instrument_files(idx).files{end+1} = fullfile(all_files(f).folder, all_files(f).name);
        instrument_files(idx).months(end+1) = month_dt;
    end
end

% For each unique instrument, read metadata from the most recent file
n_found = length(instrument_files);
for i = 1:length(instrument_files)
    instrument_id = instrument_files(i).instrument_id;
    files = instrument_files(i).files;
    available_months = unique(instrument_files(i).months);
    
    % Use first file to get metadata
    sample_file = files{1};
    
    try
        info = ncinfo(sample_file);
        
        % Initialize values
        lat = NaN;
        lon = NaN;
        station_name = '';
        instrument_type = '';
        
        % Read station_latitude, station_longitude
        var_names = {info.Variables.Name};
        if ismember('station_latitude', var_names)
            lat = double(ncread(sample_file, 'station_latitude'));
            if numel(lat) > 1, lat = lat(1); end
        end
        if ismember('station_longitude', var_names)
            lon = double(ncread(sample_file, 'station_longitude'));
            if numel(lon) > 1, lon = lon(1); end
        end
        
        % Search global attributes for site_location and instrument type
        for a = 1:length(info.Attributes)
            att_name = info.Attributes(a).Name;
            att_value = info.Attributes(a).Value;
            
            % Extract site_location
            if strcmpi(att_name, 'site_location') || strcmpi(att_name, 'location')
                if ischar(att_value) || isstring(att_value)
                    station_name = char(att_value);
                end
            end
            
            % Extract instrument type from various possible attributes
            if               strcmpi(att_name, 'instrument_type') 
                    instrument_type = char(att_value);
            end
        end
        
        % Only add if instrument type was detected
        if ~isempty(instrument_type)
            n_recognized = n_recognized + 1;
            instruments(end+1).wmo = wmo_code;
            instruments(end).instrument_id = instrument_id;
            instruments(end).instrument_type = instrument_type;
            instruments(end).available_months = available_months;
            instruments(end).lat = lat;
            instruments(end).lon = lon;
            instruments(end).station_name = station_name;
        else
            % Debug: report unrecognized instruments with raw attribute values
            fprintf('    WARNING: No instrument type found for %s instrument %s\n', wmo_code, instrument_id);
            fprintf('             File: %s\n', sample_file);
            % Print all relevant attributes for debugging
            for a = 1:length(info.Attributes)
                att_name = info.Attributes(a).Name;
                if any(strcmpi(att_name, {'source', 'instrument', 'instrument_type', 'title'}))
                    fprintf('             %s: %s\n', att_name, string(info.Attributes(a).Value));
                end
            end
        end
        
    catch ME
        fprintf('    ERROR reading %s: %s\n', sample_file, ME.message);
        continue;
    end
end

end

function [monthly_counts, instrument_types] = count_stations_per_month(all_instruments, month_dates, instrument_types)
% COUNT_STATIONS_PER_MONTH - Count unique stations with data per month per instrument type

n_months = length(month_dates);
n_types = length(instrument_types);

% Initialize count matrix: rows = months, cols = instrument types
monthly_counts = zeros(n_months, n_types);

for m = 1:n_months
    current_month = month_dates(m);
    
    for t = 1:n_types
        inst_type = instrument_types{t};
        
        % Count unique stations (wmo + instrument_id) with data this month
        count = 0;
        for i = 1:length(all_instruments)
            if strcmpi(all_instruments(i).instrument_type, inst_type)
                % Check if this instrument has data for current month
                if any(all_instruments(i).available_months == current_month)
                    count = count + 1;
                end
            end
        end
        
        monthly_counts(m, t) = count;
    end
end

end

function instrument_type = extract_instrument_type(att_str)
% EXTRACT_INSTRUMENT_TYPE - Extract instrument type from attribute string
% Simply returns the raw attribute value (trimmed)

instrument_type = strtrim(att_str);

end

function display_summary(month_dates, monthly_counts, instrument_types)
% DISPLAY_SUMMARY - Display a summary table of monthly counts

fprintf('\n========================================\n');
fprintf('         MONTHLY STATION COUNTS\n');
fprintf('========================================\n\n');

% Header
fprintf('%10s', 'Month');
for t = 1:length(instrument_types)
    fprintf('%10s', instrument_types{t});
end
fprintf('%10s\n', 'Total');
fprintf('%s\n', repmat('-', 1, 10 + 10*length(instrument_types) + 10));

% Data rows
for m = 1:length(month_dates)
    fprintf('%10s', datestr(month_dates(m), 'yyyy-mm'));
    for t = 1:length(instrument_types)
        fprintf('%10d', monthly_counts(m, t));
    end
    fprintf('%10d\n', sum(monthly_counts(m, :)));
end

% Totals
fprintf('%s\n', repmat('-', 1, 10 + 10*length(instrument_types) + 10));
fprintf('%10s', 'Max');
for t = 1:length(instrument_types)
    fprintf('%10d', max(monthly_counts(:, t)));
end
fprintf('%10d\n', max(sum(monthly_counts, 2)));

end

function plot_cumulative_bars(month_dates, monthly_counts, instrument_types, instrument_colors)
% PLOT_CUMULATIVE_BARS - Create stacked bar chart of station counts by instrument type

figure('Position', [100 100 1200 600], 'Color', 'w');

% Prepare colors (use predefined or auto-generate)
n_types = length(instrument_types);
colors = zeros(n_types, 3);
default_colormap = lines(n_types);  % Fallback colormap

for t = 1:n_types
    % Normalize instrument type for field name lookup
    field_name = regexprep(instrument_types{t}, '[^a-zA-Z0-9]', '');
    if isfield(instrument_colors, field_name)
        colors(t, :) = instrument_colors.(field_name);
    else
        colors(t, :) = default_colormap(t, :);  % Auto-assign color
    end
end

% Create stacked bar chart
h = bar(month_dates, monthly_counts, 'stacked');

% Apply colors
for t = 1:n_types
    h(t).FaceColor = colors(t, :);
    h(t).EdgeColor = 'none';
end

% Formatting
xlabel('Month', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Number of Stations', 'FontSize', 12, 'FontWeight', 'bold');
title('Stations with Valid Data by Month and Instrument Type', ...
      'FontSize', 14, 'FontWeight', 'bold');

% Legend
legend(instrument_types, 'Location', 'northwest', 'FontSize', 11);

% Grid
grid on;
box on;

% X-axis formatting
ax = gca;
ax.FontSize = 11;

% Adjust x-axis for datetime
if length(month_dates) > 12
    % For long time series, show fewer tick labels
    ax.XTick = month_dates(1:3:end);
end
ax.XTickLabelRotation = 45;



% Adjust y-axis to make room for labels
% ylim([0 max(total_counts) * 1.1]);

% Save figure
saveas(gcf, 'station_counts_by_month.png');
fprintf('\nFigure saved to: station_counts_by_month.png\n');

% Create cumulative growth chart
figure('Position', [150 150 1200 600], 'Color', 'w');

% Calculate cumulative unique stations over time
cumulative_counts = zeros(length(month_dates), n_types);
seen_stations = cell(1, n_types);  % Track unique stations per type
for t = 1:n_types
    seen_stations{t} = {};
end

% This requires access to all_instruments - need to compute differently
% For now, show cumulative sum (stations that have ever had data up to that month)
for m = 1:length(month_dates)
    cumulative_counts(m, :) = max(monthly_counts(1:m, :), [], 1);
end

% Alternative: just show the running maximum
h2 = area(month_dates, cumulative_counts);

% Apply colors
for t = 1:n_types
    h2(t).FaceColor = colors(t, :);
    h2(t).EdgeColor = 'w';
    h2(t).FaceAlpha = 0.8;
end

% Formatting
xlabel('Month', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Maximum Number of Stations', 'FontSize', 12, 'FontWeight', 'bold');
title('Peak Station Coverage Over Time by Instrument Type', ...
      'FontSize', 14, 'FontWeight', 'bold');

legend(instrument_types, 'Location', 'northwest', 'FontSize', 11);
grid on;
box on;

ax = gca;
ax.FontSize = 11;
ax.XTickLabelRotation = 45;

% Save figure
saveas(gcf, 'station_counts_cumulative.png');
fprintf('Figure saved to: station_counts_cumulative.png\n');

end
