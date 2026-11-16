% ===== Provenance (paper code archive) =====
% Paper Figure 18: climatology of first-layer cloud base height across the network (section
% 'Diurnal and Monthly Climatology - CBH', lines ~73-159 of the original file), from the E-Profile L2
% cloud data read by the local function readEProfileData() into eprofile-data.mat / eprofile-cloud-data-*.mat.
% CAVEAT (2026-10-20): the station set of the paper figure (latitude ticks 76.5, 64.9, 61.9 ... 28.5) is not
% reproduced by any MAT file now on disk; regenerate from eprofile-cloud-data-2016-2025.mat (2026-02-16).
% Source: MATLAB/ALC repository, untracked there; copied here unchanged on 2026-10-20.
% ==========================================
%% E-PROFILE DATA READER AND VISUALIZATION
% Fast parallel reading of E-Profile netCDF files
% Data structure (stationData cell array):
%   Column 1: Station ID (WIGOS_Identifier)
%   Column 2: Site Name
%   Column 3: Latitude
%   Column 4: Longitude
%   Column 5: Time [nObs x 1] datetime
%   Column 6: CBH [nObs x 3] cloud base height (m)
%   Column 7: Cloud Amount [nObs x 1] (octa)
%   Column 8: Vertical Visibility [nObs x 1] (m)
%   Column 9: Instrument Type (String)
%
% Author: Max
% Date: January 2026

clear all; close all; clc;

%% CONFIGURATION
CONFIG.folder_loc = 'A:\E-PROFILE_L2_hourly_monthly_all\';
CONFIG.outputDir = 'C:\DATA\Projects\E-PROFILE_Paper';
CONFIG.output_file =  'eprofile-cloud-data-2016-2025.mat';
CONFIG.reload_data = false; % 
CONFIG.max_files = inf; 

DATE_MIN = datetime(2016,1,1);
DATE_MAX = datetime(2025,12,31);

%% LOAD OR READ DATA
if exist(CONFIG.output_file, 'file') && ~CONFIG.reload_data
    fprintf('Loading existing data from %s...\n', CONFIG.output_file);
    load(CONFIG.output_file);
    % Check for new columns
    if size(stationData, 2) < 9
        warning('Loaded data is missing Instrument Type. Please set CONFIG.reload_data = true.');
    end
    fprintf('Loaded %d stations.\n', size(stationData, 1));
else
    [stationData, uniqueStations] = readEProfileData(CONFIG);
end

%% FILTER DATA BY DATE RANGE
fprintf('\nFiltering data to range: %s to %s\n', string(DATE_MIN), string(DATE_MAX));
for s = 1:size(stationData, 1)
    time = stationData{s, 5};
    if ~isempty(time)
        mask = time >= DATE_MIN & time <= DATE_MAX;
        stationData{s, 5} = time(mask);
        stationData{s, 6} = stationData{s, 6}(mask, :);
        stationData{s, 7} = stationData{s, 7}(mask);
        stationData{s, 8} = stationData{s, 8}(mask);
    end
end
fprintf('Data filtering complete.\n');

%% MULTI-STATION PLOTS 
nStations = size(stationData, 1);

%% PLOT: Station Map
lats = cell2mat(stationData(:,3));
lons = cell2mat(stationData(:,4));
nObs = cellfun(@length, stationData(:,5));

figure('Position', [100 100 1000 700], 'Name', 'Station Map');
gx = geoaxes;
geoscatter(gx, lats, lons, 25, nObs, 'filled');
c = colorbar;
ylabel(c, 'Number of observations');
colormap(jet);
geobasemap(gx, 'topographic');
title('E-Profile Stations');

%% PLOT: Diurnal and Monthly Climatology - CBH (All Stations)
hourlyData_cbh = NaN(nStations, 24);
monthlyData_cbh = NaN(nStations, 12);
siteNames_all = cell(nStations, 1);

for s = 1:nStations
    time_s = stationData{s, 5};
    siteNames_all{s} = stationData{s, 2};
    lon_s = stationData{s, 4};
    data_s = stationData{s, 6}(:,1);
    valid_s = data_s > 0 & ~isnan(data_s);
    
    if sum(valid_s) < 100, continue; end
    
    time_s = time_s(valid_s);
    data_s = data_s(valid_s);
    
    utc_hour = hour(time_s) + minute(time_s)/60;
    solar_hour = mod(utc_hour + lon_s/15, 24);
    solar_hour_bin = floor(solar_hour);
    
    for hr = 0:23
        mask_hr = solar_hour_bin == hr;
        if sum(mask_hr) > 10
            hourlyData_cbh(s, hr+1) = mean(data_s(mask_hr));
        end
    end
    
    for mo = 1:12
        mask_mo = month(time_s) == mo;
        if sum(mask_mo) > 10
            monthlyData_cbh(s, mo) = mean(data_s(mask_mo));
        end
    end
end

lats_sort = cell2mat(stationData(:,3));
[lats_sorted, sortIdx] = sort(lats_sort, 'descend');

hourlyData_cbh = hourlyData_cbh(sortIdx, :);
monthlyData_cbh = monthlyData_cbh(sortIdx, :);
siteNames_sorted = siteNames_all(sortIdx);

yTickStep = max(1, round(nStations/15));
yTickIdx = 1:yTickStep:nStations;

figure('Position', [100 100 1000 700], 'Name', 'Cloud Base Height Climatology');

subplot(1,2,1);
hourlyData_cbh_ext = [hourlyData_cbh, hourlyData_cbh(:,1)];
hourlyData_cbh_ext = [hourlyData_cbh_ext; hourlyData_cbh_ext(end,:)];
pcolor(0:24, 1:nStations+1, hourlyData_cbh_ext);
shading flat;
colorbar;
colormap(gca, parula);
clim([0 3000]);
xlabel('Hour (Solar Time)', 'FontSize', 11);
ylabel('Latitude (°N)', 'FontSize', 11);
title('(a) Diurnal Cycle', 'FontSize', 12, 'FontWeight', 'normal');
yticks(yTickIdx + 0.5);
yticklabels(compose('%.1f', lats_sorted(yTickIdx)));
xticks(0:6:24);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlim([0 24]);
ylim([1 nStations+1]);

subplot(1,2,2);
monthlyData_cbh_ext = [monthlyData_cbh, NaN(nStations,1)];
monthlyData_cbh_ext = [monthlyData_cbh_ext; monthlyData_cbh_ext(end,:)];
pcolor(1:13, 1:nStations+1, monthlyData_cbh_ext);
shading flat;
cb = colorbar;
cb.Label.String = 'CBH (m)';
cb.Label.FontSize = 11;
colormap(gca, parula);
clim([0 3000]);
xlabel('Month', 'FontSize', 11);
title('(b) Seasonal Cycle', 'FontSize', 12, 'FontWeight', 'normal');
yticks(yTickIdx + 0.5);
yticklabels({});
xticks((1:12) + 0.5);
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'});
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlim([1 13]);
ylim([1 nStations+1]);

sgtitle('Cloud Base Height Climatology (Layer 1)', 'FontSize', 13, 'FontWeight', 'bold');

%% PLOT: Diurnal and Monthly Climatology - Cloud Amount (All Stations)
hourlyData_cloud = NaN(nStations, 24);
monthlyData_cloud = NaN(nStations, 12);

for s = 1:nStations
    time_s = stationData{s, 5};
    lon_s = stationData{s, 4};
    data_s = stationData{s, 7};
    valid_s = ~isnan(data_s);
    
    if sum(valid_s) < 100, continue; end
    
    time_s = time_s(valid_s);
    data_s = data_s(valid_s);
    
    utc_hour = hour(time_s) + minute(time_s)/60;
    solar_hour = mod(utc_hour + lon_s/15, 24);
    solar_hour_bin = floor(solar_hour);
    
    for hr = 0:23
        mask_hr = solar_hour_bin == hr;
        if sum(mask_hr) > 10
            hourlyData_cloud(s, hr+1) = mean(data_s(mask_hr));
        end
    end
    
    for mo = 1:12
        mask_mo = month(time_s) == mo;
        if sum(mask_mo) > 10
            monthlyData_cloud(s, mo) = mean(data_s(mask_mo));
        end
    end
end

hourlyData_cloud = hourlyData_cloud(sortIdx, :);
monthlyData_cloud = monthlyData_cloud(sortIdx, :);

figure('Position', [100 100 1000 700], 'Name', 'Cloud Amount Climatology');

subplot(1,2,1);
hourlyData_cloud_ext = [hourlyData_cloud, hourlyData_cloud(:,1)];
hourlyData_cloud_ext = [hourlyData_cloud_ext; hourlyData_cloud_ext(end,:)];
pcolor(0:24, 1:nStations+1, hourlyData_cloud_ext);
shading flat;
colorbar;
colormap(gca, flipud(gray));
clim([0 8]);
xlabel('Hour (Solar Time)', 'FontSize', 11);
ylabel('Latitude (°N)', 'FontSize', 11);
title('(a) Diurnal Cycle', 'FontSize', 12, 'FontWeight', 'normal');
yticks(yTickIdx + 0.5);
yticklabels(compose('%.1f', lats_sorted(yTickIdx)));
xticks(0:6:24);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlim([0 24]);
ylim([1 nStations+1]);

subplot(1,2,2);
monthlyData_cloud_ext = [monthlyData_cloud, NaN(nStations,1)];
monthlyData_cloud_ext = [monthlyData_cloud_ext; monthlyData_cloud_ext(end,:)];
pcolor(1:13, 1:nStations+1, monthlyData_cloud_ext);
shading flat;
cb = colorbar;
cb.Label.String = 'Cloud Amount (octa)';
cb.Label.FontSize = 11;
colormap(gca, flipud(gray));
clim([0 8]);
xlabel('Month', 'FontSize', 11);
title('(b) Seasonal Cycle', 'FontSize', 12, 'FontWeight', 'normal');
yticks(yTickIdx + 0.5);
yticklabels({});
xticks((1:12) + 0.5);
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'});
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlim([1 13]);
ylim([1 nStations+1]);

sgtitle('Cloud Amount Climatology', 'FontSize', 13, 'FontWeight', 'bold');

%% PLOT: Monthly Time Series Heatmaps
monthBins = dateshift(DATE_MIN, 'start', 'month') : calmonths(1) : dateshift(DATE_MAX, 'end', 'month');
nMonths = length(monthBins) - 1;

% Cloud Amount Time Series
monthlyData_ts = NaN(nStations, nMonths);
validMask_ts = false(nStations, 1);

for s = 1:nStations
    time_s = stationData{s, 5};
    data_s = stationData{s, 7};
    valid_s = ~isnan(data_s) & data_s >= 0;
    if sum(valid_s) < 100, continue; end
    validMask_ts(s) = true;
    time_s = time_s(valid_s);
    data_s = data_s(valid_s);
    for m = 1:nMonths
        mask_m = time_s >= monthBins(m) & time_s < monthBins(m+1);
        if sum(mask_m) > 10
            monthlyData_ts(s, m) = mean(data_s(mask_m));
        end
    end
end

monthlyData_ts = monthlyData_ts(validMask_ts, :);
lats_ts = cell2mat(stationData(validMask_ts, 3));
[lats_ts_sorted, sortIdx_ts] = sort(lats_ts, 'descend');
monthlyData_ts = monthlyData_ts(sortIdx_ts, :);
nValidStations = size(monthlyData_ts, 1);
yTickStep_ts = max(1, round(nValidStations/15));
yTickIdx_ts = 1:yTickStep_ts:nValidStations;

figure('Position', [100 100 1400 600], 'Name', 'Monthly Mean Cloud Cover');
monthlyData_ts_ext = [monthlyData_ts, NaN(nValidStations,1)];
monthlyData_ts_ext = [monthlyData_ts_ext; monthlyData_ts_ext(end,:)];
pcolor(monthBins, 1:nValidStations+1, monthlyData_ts_ext);
shading flat;
cb = colorbar;
cb.Label.String = 'Cloud Amount (octa)';
cb.Label.FontSize = 11;
colormap(flipud(parula));
clim([0 8]);
xlabel('Time', 'FontSize', 11);
ylabel('Latitude (°N)', 'FontSize', 11);
title('Monthly Mean Cloud Cover', 'FontSize', 13);
yticks(yTickIdx_ts + 0.5);
yticklabels(compose('%.1f', lats_ts_sorted(yTickIdx_ts)));
set(gca, 'YDir', 'reverse', 'FontSize', 10);
ylim([1 nValidStations+1]);

% CBH Time Series
monthlyData_cbh_ts = NaN(nStations, nMonths);
validMask_cbh_ts = false(nStations, 1);

for s = 1:nStations
    time_s = stationData{s, 5};
    data_s = stationData{s, 6}(:,1);
    valid_s = data_s > 0 & ~isnan(data_s);
    if sum(valid_s) < 100, continue; end
    validMask_cbh_ts(s) = true;
    time_s = time_s(valid_s);
    data_s = data_s(valid_s);
    for m = 1:nMonths
        mask_m = time_s >= monthBins(m) & time_s < monthBins(m+1);
        if sum(mask_m) > 10
            monthlyData_cbh_ts(s, m) = mean(data_s(mask_m));
        end
    end
end

monthlyData_cbh_ts = monthlyData_cbh_ts(validMask_cbh_ts, :);
lats_cbh_ts = cell2mat(stationData(validMask_cbh_ts, 3));
[lats_cbh_ts_sorted, sortIdx_cbh_ts] = sort(lats_cbh_ts, 'descend');
monthlyData_cbh_ts = monthlyData_cbh_ts(sortIdx_cbh_ts, :);
nValidStations_cbh = size(monthlyData_cbh_ts, 1);
yTickStep_cbh = max(1, round(nValidStations_cbh/15));
yTickIdx_cbh = 1:yTickStep_cbh:nValidStations_cbh;

figure('Position', [100 100 1400 600], 'Name', 'Monthly Mean Cloud Base Height');
monthlyData_cbh_ts_ext = [monthlyData_cbh_ts, NaN(nValidStations_cbh,1)];
monthlyData_cbh_ts_ext = [monthlyData_cbh_ts_ext; monthlyData_cbh_ts_ext(end,:)];
pcolor(monthBins, 1:nValidStations_cbh+1, monthlyData_cbh_ts_ext);
shading flat;
cb = colorbar;
cb.Label.String = 'CBH (m)';
cb.Label.FontSize = 11;
colormap(parula);
clim([0 3000]);
xlabel('Time', 'FontSize', 11);
ylabel('Latitude (°N)', 'FontSize', 11);
title('Monthly Mean Cloud Base Height (Layer 1)', 'FontSize', 13);
yticks(yTickIdx_cbh + 0.5);
yticklabels(compose('%.1f', lats_cbh_ts_sorted(yTickIdx_cbh)));
set(gca, 'YDir', 'reverse', 'FontSize', 10);
ylim([1 nValidStations_cbh+1]);

%% PLOT: Data Availability
monthlyCount = zeros(nStations, nMonths);
fprintf('Calculating monthly observations...\n');
tic
for s = 1:nStations
    monthlyCount(s,:) = histcounts(stationData{s, 5}, monthBins);
end
toc

lats_avail = cell2mat(stationData(:,3));
[lats_avail_sorted, sortIdx_avail] = sort(lats_avail, 'descend');
monthlyCount = monthlyCount(sortIdx_avail, :);

% EXTRACT METADATA
avail_ids = stationData(sortIdx_avail, 1);    % WIGOS_Identifier
avail_names = stationData(sortIdx_avail, 2);  % Name
avail_types = stationData(sortIdx_avail, 9);  % Instrument Type

% --- CONSTRUCT LABELS ---
% Format: Name, WIGOS_ID [Identifier] Type
fullLabels = cell(length(avail_ids), 1);
for k = 1:length(avail_ids)
    % Parse WIGOS and Identifier
    parts = split(avail_ids{k}, '_');
    wigos = parts{1};
    if length(parts) > 1
        ident = parts{2}; 
    else
        ident = '-'; 
    end
    
    name = avail_names{k};
    
    instType = avail_types{k};
    if isempty(instType), instType = '?'; end
    
    % Combine: "Payerne, 0-20000-0-06610 [A] CHM15k"
    fullLabels{k} = sprintf('%s, %s [%s] %s', replace(name,'_',' '), wigos, ident, instType);
end

maxPerPanel = 60; 
nTotalPanels = ceil(nStations / maxPerPanel);
dataLog = log10(monthlyCount + 1);
cLimits = [0, max(dataLog(:))];

for p = 1:nTotalPanels
    idxStart = (p - 1) * maxPerPanel + 1;
    idxEnd = min(p * maxPerPanel, nStations);
    panelIdx = idxStart:idxEnd;
    nPanelStations = length(panelIdx);
    
    % Widen figure for long labels
    fig_h = figure('Position', [50 50 1400 950], ... 
           'Name', sprintf('Data Availability - Stations %d-%d', idxStart, idxEnd));
    
    % Adjust axes: More space on left [left bottom width height]
    ax_p = axes('Position', [0.45 0.08 0.35 0.85]); 

    dataLog_panel = dataLog(panelIdx, :);
    dataLog_panel_ext = [dataLog_panel, NaN(nPanelStations,1)];
    dataLog_panel_ext = [dataLog_panel_ext; dataLog_panel_ext(end,:)];
    
    h_pcolor = pcolor(monthBins, 1:nPanelStations+1, dataLog_panel_ext);
    shading flat;
    ax_p.Tag = sprintf('panel_%d', p);
    
    % Store data for datatip
    panelData = struct();
    panelData.siteNames = avail_names(panelIdx);
    panelData.monthBins = monthBins;
    panelData.stationIDs = avail_ids(panelIdx);
    panelData.instTypes = avail_types(panelIdx); 
    panelData.obsCount = monthlyCount(panelIdx, :);
    panelData.idxStart = idxStart;
    
    clim(cLimits);
    colormap(sky);
    
    % Left Y-axis: Full Labels
    yTickIdx_p = 1:nPanelStations;
    yticks(yTickIdx_p + 0.5);
    yticklabels(fullLabels(panelIdx));
    
    set(ax_p, 'YDir', 'reverse', 'FontSize', 8); 
    ylim([1 nPanelStations+1]);
    ylabel('Stations', 'FontSize', 11);
    xlabel('Time', 'FontSize', 11);
    
    % Right Y-axis: Latitude
    ax_right = axes('Position', ax_p.Position, ...
        'YAxisLocation', 'right', ...
        'Color', 'none', ...
        'XTick', [], ...
        'YDir', 'reverse', ...
        'FontSize', 7);
    ax_right.YLim = [1 nPanelStations+1];
    ax_right.YTick = yTickIdx_p + 0.5;
    ax_right.YTickLabel = compose('%.2f°', lats_avail_sorted(panelIdx));
    ylabel(ax_right, 'Latitude (°N)', 'FontSize', 11);
    
    linkaxes([ax_p, ax_right], 'y');
    
    cb = colorbar(ax_p);
    cb.Position = [0.85 0.08 0.02 0.85];
    cb.Label.String = 'log_{10}(Observations)';
    
    title(ax_p, sprintf('Data Availability - Stations %d to %d', idxStart, idxEnd), 'FontSize', 12);
    setappdata(fig_h, 'panelData', {panelData});
    dcm = datacursormode(fig_h);
    set(dcm, 'UpdateFcn', @(src,evt) availabilityDataTip(src, evt, fig_h));
end

fprintf('Created %d Data Availability figure(s) for %d stations\n', nTotalPanels, nStations);

%% SAVE ALL FIGURES
if ~exist(CONFIG.outputDir, 'dir')
    mkdir(CONFIG.outputDir);
end
fprintf('\nSaving figures to %s...\n', CONFIG.outputDir);
figs = findall(0, 'Type', 'figure');
for i = 1:length(figs)
    fig = figs(i);
    figName = get(fig, 'Name');
    if isempty(figName), figName = sprintf('Figure_%d', fig.Number); end
    figName = regexprep(figName, '[\\/:*?"<>|]', '_');
    figName = strrep(figName, ' ', '_');
    pngFile = fullfile(CONFIG.outputDir, [figName '.png']);
    fprintf('  Saving: %s\n', figName);
    exportgraphics(fig, pngFile, 'Resolution', 300);
end
fprintf('Done! Saved %d figures.\n', length(figs));

%%  INTERACTIVE STATION BROWSER
siteNames_browser = stationData(:, 2);
stationIDs_browser = stationData(:, 1);

labels_browser = cell(nStations, 1);
for i = 1:nStations
    labels_browser{i} = sprintf('%s (%s)', siteNames_browser{i}, stationIDs_browser{i});
end
[labels_sorted, sortIdx_browser] = sort(labels_browser);

fig_browser = figure('Position', [50 50 1800 950], ...
    'Name', 'Interactive Station Browser', 'Color', 'w', 'NumberTitle', 'off');

uicontrol('Style', 'text', 'String', 'Select Station:', ...
    'Position', [20 915 100 20], 'FontSize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');

dd = uicontrol('Style', 'popupmenu', 'String', labels_sorted, ...
    'Position', [120 910 350 30], 'FontSize', 10, 'Callback', @updateAllPlots);

infoText = uicontrol('Style', 'text', 'String', '', ...
    'Position', [500 905 800 30], 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', 'w');

t_browser = tiledlayout(fig_browser, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
t_browser.OuterPosition = [0 0 1 0.95];

ax_cbh_ts = nexttile(t_browser, 1);       
ax_ca_ts = nexttile(t_browser, 2);         
ax_density = nexttile(t_browser, 3);       
ax_cbh_hist = nexttile(t_browser, 4);      
ax_ca_hist = nexttile(t_browser, 5);       
ax_diurnal = nexttile(t_browser, 6);       
ax_monthly_cbh = nexttile(t_browser, 7);   
ax_monthly_ca = nexttile(t_browser, 8);    
ax_combined = nexttile(t_browser, 9);      

setappdata(fig_browser, 'stationData', stationData);
setappdata(fig_browser, 'sortIdx_browser', sortIdx_browser);
setappdata(fig_browser, 'DATE_MIN', DATE_MIN);
setappdata(fig_browser, 'DATE_MAX', DATE_MAX);
setappdata(fig_browser, 'axes', struct('cbh_ts', ax_cbh_ts, 'ca_ts', ax_ca_ts, ...
    'density', ax_density, 'cbh_hist', ax_cbh_hist, 'ca_hist', ax_ca_hist, ...
    'diurnal', ax_diurnal, 'monthly_cbh', ax_monthly_cbh, ...
    'monthly_ca', ax_monthly_ca, 'combined', ax_combined));
setappdata(fig_browser, 'infoText', infoText);

updateAllPlots(dd, []);

function updateAllPlots(src, ~)
    fig = ancestor(src, 'figure');
    stationData_local = getappdata(fig, 'stationData');
    DATE_MIN_local = getappdata(fig, 'DATE_MIN');
    DATE_MAX_local = getappdata(fig, 'DATE_MAX');
    ax = getappdata(fig, 'axes');
    infoTextHandle = getappdata(fig, 'infoText');
    sortIdx_browser = getappdata(fig, 'sortIdx_browser');
    
    idx_s = sortIdx_browser(src.Value);
    
    stationID_s = stationData_local{idx_s, 1};
    lat_s = stationData_local{idx_s, 3};
    lon_s = stationData_local{idx_s, 4};
    time_s = stationData_local{idx_s, 5};
    cbh_s = stationData_local{idx_s, 6};
    ca_s = stationData_local{idx_s, 7};
    
    if size(stationData_local, 2) >= 9
        instType = stationData_local{idx_s, 9};
    else
        instType = 'N/A';
    end

    infoTextHandle.String = sprintf('ID: %s | Type: %s | Lat: %.3f | Lon: %.3f | Obs: %d', ...
        stationID_s, instType, lat_s, lon_s, length(time_s));
    
    % PLOT 1: CBH Time Series 
    axes(ax.cbh_ts); cla; hold on;
    colors_cbh = {[0 0.45 0.74], [0.85 0.33 0.1], [0.93 0.69 0.13]};
    layerNames = {'L1', 'L2', 'L3'};
    for layer = 1:3
        cbh_layer = cbh_s(:, layer);
        valid_layer = cbh_layer > 0 & ~isnan(cbh_layer);
        if sum(valid_layer) > 0
            plot(time_s(valid_layer), cbh_layer(valid_layer), '.', 'MarkerSize', 1, 'Color', colors_cbh{layer});
        end
    end
    hold off; ylabel('CBH (m)'); title('Cloud Base Height Time Series');
    legend(layerNames, 'Location', 'northeast', 'FontSize', 8);
    ylim([0 12000]); xlim([DATE_MIN_local, DATE_MAX_local]); grid on; box on;
    
    % PLOT 2: Cloud Amount Time Series
    axes(ax.ca_ts); cla; hold on;
    valid_ca = ~isnan(ca_s);
    if sum(valid_ca) > 0
        plot(time_s(valid_ca), ca_s(valid_ca), '.', 'MarkerSize', 1);
        TT = timetable(time_s(valid_ca), ca_s(valid_ca), 'VariableNames', {'CA'});
        TT_monthly = retime(TT, 'monthly', @median);
        plot(TT_monthly.Time, TT_monthly.CA, '-', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
    end
    hold off; ylabel('Cloud Amount (octa)'); title('Cloud Amount Time Series');
    ylim([0 9]); xlim([DATE_MIN_local, DATE_MAX_local]); grid on; box on;
    legend({'Raw', 'Monthly median'}, 'Location', 'best', 'FontSize', 8);
    linkaxes([ax.cbh_ts, ax.ca_ts], 'x');
    
    % PLOT 3: CBH Density Heatmap
    axes(ax.density); cla;
    cbh1 = cbh_s(:,1); mask_heat = cbh1 > 0 & ~isnan(cbh1);
    time_heat = time_s(mask_heat); cbh1_heat = cbh1(mask_heat);
    if ~isempty(time_heat)
        timeBins = DATE_MIN_local:days(7):DATE_MAX_local; 
        heightBins = 0:100:10000;
        [~, ~, tBin] = histcounts(time_heat, timeBins);
        [~, ~, hBin] = histcounts(cbh1_heat, heightBins);
        valid_heat = tBin > 0 & hBin > 0;
        if sum(valid_heat) > 0
            heatmapData = accumarray([hBin(valid_heat), tBin(valid_heat)], 1, [length(heightBins)-1, length(timeBins)-1]);
            heatmapData_ext = [heatmapData, zeros(size(heatmapData,1),1)];
            heatmapData_ext = [heatmapData_ext; zeros(1, size(heatmapData_ext,2))];
            pcolor(timeBins, heightBins, heatmapData_ext); shading flat;
            colorbar; colormap(ax.density, sky); clim([0 prctile(heatmapData(:), 99)]);
        end
    end
    xlabel('Date'); ylabel('CBH (m)'); title('CBH Density (Layer 1)'); ylim([0 6000]);
    
    % PLOT 4: CBH Histogram
    axes(ax.cbh_hist); cla; hold on;
    colors_hist = {[0 0.45 0.74], [0.85 0.33 0.1], [0.93 0.69 0.13]};
    for layer = 1:3
        data_hist = cbh_s(:, layer);
        data_hist = data_hist(data_hist > 0 & ~isnan(data_hist));
        if ~isempty(data_hist), histogram(data_hist, 0:200:10000, 'FaceColor', colors_hist{layer}, 'FaceAlpha', 0.5); end
    end
    hold off; xlabel('CBH (m)'); ylabel('Count'); xlim([0 10000]);
    legend(layerNames, 'Location', 'northeast', 'FontSize', 8); title('CBH Distribution'); grid on;
    
    % PLOT 5: Cloud Amount Histogram
    axes(ax.ca_hist); cla;
    valid_ca_hist = ca_s(~isnan(ca_s) & ca_s >= 0);
    if ~isempty(valid_ca_hist), histogram(valid_ca_hist, -0.5:1:8.5, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'w'); end
    xlabel('Cloud Amount (octa)'); ylabel('Count'); xlim([-0.5 8.5]); xticks(0:8); title('Cloud Amount Distribution'); grid on;
    
    % PLOT 6: Diurnal Cycle
    axes(ax.diurnal); cla;
    cbh1_diurnal = cbh_s(:,1); valid_diurnal = cbh1_diurnal > 0 & ~isnan(cbh1_diurnal);
    if sum(valid_diurnal) > 100
        time_diurnal = time_s(valid_diurnal); cbh1_diurnal = cbh1_diurnal(valid_diurnal);
        utc_hour_diurnal = hour(time_diurnal) + minute(time_diurnal)/60;
        solar_hour_diurnal = mod(utc_hour_diurnal + lon_s/15, 24);
        solar_hour_bin = floor(solar_hour_diurnal);
        hourlyMean = accumarray(solar_hour_bin+1, cbh1_diurnal, [24 1], @mean);
        hourlyStd = accumarray(solar_hour_bin+1, cbh1_diurnal, [24 1], @std);
        hourlyCount = accumarray(solar_hour_bin+1, cbh1_diurnal, [24 1], @numel);
        hourlyCount(hourlyCount == 0) = NaN; hourlyMean(hourlyCount < 10) = NaN;
        errorbar(0:23, hourlyMean, hourlyStd./sqrt(hourlyCount), '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0 0.45 0.74]);
    end
    xlabel('Hour (Solar Time)'); ylabel('Mean CBH L1 (m)'); title('Diurnal Cycle'); grid on; xlim([-0.5 23.5]); xticks(0:6:23);
    
    % PLOT 7: Monthly CBH
    axes(ax.monthly_cbh); cla;
    cbh1_monthly = cbh_s(:,1); valid_monthly = cbh1_monthly > 0 & ~isnan(cbh1_monthly);
    if sum(valid_monthly) > 100
        time_monthly = time_s(valid_monthly); cbh1_monthly = cbh1_monthly(valid_monthly);
        monthlyMean_cbh = zeros(12, 1); monthlyStd_cbh = zeros(12, 1); monthlyCount_cbh = zeros(12, 1);
        for mo = 1:12
            mask_mo = month(time_monthly) == mo;
            if sum(mask_mo) > 10
                monthlyMean_cbh(mo) = mean(cbh1_monthly(mask_mo));
                monthlyStd_cbh(mo) = std(cbh1_monthly(mask_mo));
                monthlyCount_cbh(mo) = sum(mask_mo);
            else, monthlyMean_cbh(mo) = NaN; end
        end
        monthlyCount_cbh(monthlyCount_cbh == 0) = NaN;
        bar(1:12, monthlyMean_cbh, 'FaceColor', [0 0.45 0.74], 'EdgeColor', 'none'); hold on;
        errorbar(1:12, monthlyMean_cbh, monthlyStd_cbh./sqrt(monthlyCount_cbh), 'k.', 'LineWidth', 1); hold off;
    end
    xlabel('Month'); ylabel('Mean CBH L1 (m)'); title('Monthly CBH Climatology'); xlim([0.5 12.5]); xticks(1:12);
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'}); grid on;
    
    % PLOT 8: Monthly Cloud Amount
    axes(ax.monthly_ca); cla;
    valid_ca_monthly = ~isnan(ca_s) & ca_s >= 0;
    if sum(valid_ca_monthly) > 100
        time_ca_monthly = time_s(valid_ca_monthly); ca_monthly = ca_s(valid_ca_monthly);
        monthlyMean_ca = zeros(12, 1); monthlyStd_ca = zeros(12, 1); monthlyCount_ca = zeros(12, 1);
        for mo = 1:12
            mask_mo = month(time_ca_monthly) == mo;
            if sum(mask_mo) > 10
                monthlyMean_ca(mo) = mean(ca_monthly(mask_mo));
                monthlyStd_ca(mo) = std(ca_monthly(mask_mo));
                monthlyCount_ca(mo) = sum(mask_mo);
            else, monthlyMean_ca(mo) = NaN; end
        end
        monthlyCount_ca(monthlyCount_ca == 0) = NaN;
        bar(1:12, monthlyMean_ca, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); hold on;
        errorbar(1:12, monthlyMean_ca, monthlyStd_ca./sqrt(monthlyCount_ca), 'k.', 'LineWidth', 1); hold off;
    end
    xlabel('Month'); ylabel('Mean Cloud Amount (octa)'); title('Monthly Cloud Amount Climatology'); xlim([0.5 12.5]); ylim([0 8]);
    xticks(1:12); xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'}); grid on;
    
    % PLOT 9: Combined Heatmap
    axes(ax.combined); cla;
    cbh1_combined = cbh_s(:,1); valid_combined = cbh1_combined > 0 & ~isnan(cbh1_combined);
    if sum(valid_combined) > 100
        time_combined = time_s(valid_combined); cbh1_combined = cbh1_combined(valid_combined);
        utc_hour_combined = hour(time_combined) + minute(time_combined)/60;
        solar_hour_combined = mod(utc_hour_combined + lon_s/15, 24); solar_hour_bin_combined = floor(solar_hour_combined);
        month_combined = month(time_combined);
        heatmap_combined = NaN(24, 12);
        for hr = 0:23
            for mo = 1:12
                mask_hm = solar_hour_bin_combined == hr & month_combined == mo;
                if sum(mask_hm) > 5, heatmap_combined(hr+1, mo) = mean(cbh1_combined(mask_hm)); end
            end
        end
        heatmap_ext = [heatmap_combined, NaN(24,1)]; heatmap_ext = [heatmap_ext; heatmap_ext(end,:)];
        pcolor(1:13, 0:24, heatmap_ext); shading flat; cb = colorbar; cb.Label.String = 'CBH (m)';
        colormap(ax.combined, parula); clim([0 prctile(heatmap_combined(:), 95)]);
    end
    xlabel('Month'); ylabel('Hour (Solar Time)'); title('Diurnal-Monthly CBH Pattern');
    xlim([1 13]); ylim([0 24]); xticks((1:12) + 0.5); xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'}); yticks(0:6:24);
    drawnow;
end

%% FUNCTION: READ E-PROFILE DATA
function [stationData, uniqueStations] = readEProfileData(CONFIG)
formatTime = @(s) char(duration(0,0,s, 'Format', 'hh:mm:ss'));
if isempty(gcp('nocreate')), parpool('local'); end
pool = gcp('nocreate'); fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
fprintf('Scanning for NetCDF files...\n');
files = dir([CONFIG.folder_loc '**/*.nc']); files = files(~[files.isdir]);
filePaths = fullfile({files.folder}, {files.name})'; nFiles = length(filePaths);
fprintf('Found %d NetCDF files\n', nFiles);
if isfield(CONFIG, 'max_files') && CONFIG.max_files < nFiles
    fprintf('Limiting to %d files for testing\n', CONFIG.max_files);
    filePaths = filePaths(1:CONFIG.max_files); nFiles = CONFIG.max_files;
end
fprintf('\nReading all files in parallel...\n');
tic;
% Added inst_type_all
[stationID_all, siteName_all, time_all, cbh_all, cloud_amount_all, vert_vis_all, lat_all, lon_all, inst_type_all] = deal(cell(nFiles, 1));
progressQueue = parallel.pool.DataQueue; progressCount = 0; lastPercent = 0; startTime = tic;
afterEach(progressQueue, @updateProgress);
    function updateProgress(~)
        progressCount = progressCount + 1; currentPercent = floor(progressCount / nFiles * 100);
        if currentPercent > lastPercent
            lastPercent = currentPercent; elapsed = toc(startTime); rate = progressCount / elapsed;
            remaining = (nFiles - progressCount) / rate; barLength = 30;
            filled = round(barLength * progressCount / nFiles);
            bar = [repmat('#', 1, filled), repmat('-', 1, barLength - filled)];
            fprintf('\r[%s] %3d%% | %d/%d files | %.1f files/s | ETA: %s', bar, currentPercent, progressCount, nFiles, rate, formatTime(remaining));
        end
    end
parfor f = 1:nFiles
    try
        [~, fname] = fileparts(filePaths{f}); tmp = split(fname, '_');
        if length(tmp) >= 3, stationID_all{f} = [tmp{2} '_' tmp{3}(1)]; else, stationID_all{f} = 'Unknown'; end
        time_raw = ncread(filePaths{f}, 'time');
        if time_raw<1e5, time_all{f} = datetime(1970,1,1) + days(time_raw); else, time_all{f} = datetime(time_raw,'ConvertFrom','datenum') ; end
        cbh = ncread(filePaths{f}, 'cloud_base_height'); cbh_all{f} = cbh';
        try cloud_amount_all{f} = ncread(filePaths{f}, 'cloud_amount'); catch, cloud_amount_all{f} = NaN(length(time_raw), 1); end
        try vert_vis_all{f} = ncread(filePaths{f}, 'vertical_visibility'); catch, vert_vis_all{f} = NaN(length(time_raw), 1); end
        try siteName_all{f} = ncreadatt(filePaths{f}, '/', 'site_location'); lat_all{f} = ncread(filePaths{f}, 'station_latitude'); lon_all{f} = ncread(filePaths{f}, 'station_longitude');
        catch, siteName_all{f} = 'Unknown'; lat_all{f} = NaN; lon_all{f} = NaN; end
        % READ INSTRUMENT TYPE
        try inst_type_all{f} = ncreadatt(filePaths{f}, '/', 'instrument_type'); catch, inst_type_all{f} = ''; end
    catch, stationID_all{f} = ''; end
    send(progressQueue, f);
end
elapsed = toc; fprintf('\nRead completed in %s (%.1f files/sec)\n', formatTime(elapsed), nFiles/elapsed);
validMask = ~cellfun(@isempty, stationID_all);
stationID_all = stationID_all(validMask); siteName_all = siteName_all(validMask);
time_all = time_all(validMask); cbh_all = cbh_all(validMask);
cloud_amount_all = cloud_amount_all(validMask); vert_vis_all = vert_vis_all(validMask);
lat_all = lat_all(validMask); lon_all = lon_all(validMask); inst_type_all = inst_type_all(validMask);
fprintf('\nAggregating by station...\n');
uniqueStations = unique(stationID_all); nStationsLocal = length(uniqueStations);
fprintf('Found %d unique stations\n\n', nStationsLocal);
stationData = cell(nStationsLocal, 9);
totalObs = 0; totalDupes = 0;
fprintf('%-15s %-20s %8s %10s %10s %15s\n', 'Station ID', 'Site Name', 'Files', 'Raw Obs', 'Dupes Rem', 'Date Range');
fprintf('%s\n', repmat('-', 1, 90));
for s = 1:nStationsLocal
    sid = uniqueStations{s}; mask = strcmp(stationID_all, sid); nFilesStation = sum(mask);
    time_cat = vertcat(time_all{mask}); cbh_cat = vertcat(cbh_all{mask});
    ca_cat = vertcat(cloud_amount_all{mask}); vv_cat = vertcat(vert_vis_all{mask});
    nRaw = length(time_cat); [time_cat, sortIdxLocal] = sort(time_cat);
    cbh_cat = cbh_cat(sortIdxLocal, :); ca_cat = ca_cat(sortIdxLocal); vv_cat = vv_cat(sortIdxLocal);
    [time_cat, uniqueIdx] = unique(time_cat);
    cbh_cat = cbh_cat(uniqueIdx, :); ca_cat = ca_cat(uniqueIdx); vv_cat = vv_cat(uniqueIdx);
    nClean = length(time_cat); nDupes = nRaw - nClean; totalDupes = totalDupes + nDupes; idxFirst = find(mask, 1);
    stationData{s, 1} = sid; stationData{s, 2} = siteName_all{idxFirst};
    stationData{s, 3} = lat_all{idxFirst}; stationData{s, 4} = lon_all{idxFirst};
    stationData{s, 5} = time_cat; stationData{s, 6} = cbh_cat;
    stationData{s, 7} = ca_cat; stationData{s, 8} = vv_cat;
    stationData{s, 9} = inst_type_all{idxFirst};
    totalObs = totalObs + nClean;
    if ~isempty(time_cat), dateRangeStr = sprintf('%s - %s', string(min(time_cat), 'yyyy-MM'), string(max(time_cat), 'yyyy-MM')); else, dateRangeStr = 'No Data'; end
    fprintf('%-15s %-20s %8d %10d %10d %15s\n', sid, extractBefore(siteName_all{idxFirst}, min(length(siteName_all{idxFirst})+1, 20)), nFilesStation, nRaw, nDupes, dateRangeStr);
end
fprintf('%s\n', repmat('-', 1, 90)); fprintf('%-15s %-20s %8d %10d %10d\n', 'TOTAL', '', nFiles, totalObs + totalDupes, totalDupes);
fprintf('\nSaving to %s...\n', CONFIG.output_file); save(CONFIG.output_file, 'stationData', 'uniqueStations', '-v7.3');
fileInfo = dir(CONFIG.output_file); fprintf('Done! File size: %.1f MB\n', fileInfo.bytes / 1e6);
end

function txt = availabilityDataTip(~, event, fig_h)
    try
        pos = event.Position;
        panelData = getappdata(fig_h, 'panelData');
        if iscell(panelData), pd = panelData{1}; else, pd = panelData; end
        if isempty(pd), txt = {sprintf('Date: %s', datestr(pos(1), 'mmm-yyyy')), sprintf('Row: %d', floor(pos(2)))}; return; end
        stationIdx = floor(pos(2));
        posTime = datetime(pos(1), 'ConvertFrom', 'datenum');
        monthIdx = discretize(posTime, pd.monthBins);
        if isempty(monthIdx) || isnan(monthIdx), monthIdx = 1; end
        if stationIdx >= 1 && stationIdx <= size(pd.obsCount, 1)
            siteName = pd.siteNames{stationIdx}; stationID = pd.stationIDs{stationIdx};
            instType = pd.instTypes{stationIdx};
            parts = split(stationID, '_'); wigosID = parts{1};
            if length(parts) > 1, identifier = parts{2}; else, identifier = ''; end
            monthIdx = min(monthIdx, size(pd.obsCount, 2)); obs = pd.obsCount(stationIdx, monthIdx);
            monthStr = datestr(pd.monthBins(monthIdx), 'mmm-yyyy');
            txt = {sprintf('Station: %s', siteName), sprintf('WIGOS ID: %s', wigosID), sprintf('Identifier: %s', identifier), sprintf('Type: %s', instType), sprintf('Date: %s', monthStr), sprintf('Obs: %d', obs)};
        else
            txt = {sprintf('Date: %s', datestr(pos(1), 'mmm-yyyy')), sprintf('Row: %d', stationIdx)};
        end
    catch ME, txt = {sprintf('X: %.1f', event.Position(1)), sprintf('Y: %.1f', event.Position(2)), sprintf('Error: %s', ME.message)}; end
end