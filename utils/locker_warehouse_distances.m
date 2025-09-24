%% locker_warehouse_distances.m
% Usage:
%   results = locker_warehouse_distances('lockers.xlsx','warehouses.xlsx');
%   % results.D is n_lockers x n_warehouses distance matrix (km)
%
% The function will save "distance_matrix.csv" and show a plot.

function results = locker_warehouse_distances(lockerFile, warehouseFile)
    if nargin<2
        error('Provide locker and warehouse file names (Excel or text CSV).');
    end

    % --------- Read tables (Excel or text CSV/TSV with common delimiters) ----------
    lockers = read_table_flexible(lockerFile);
    warehouses = read_table_flexible(warehouseFile);

    % Normalize and get lat/lon columns (case-insensitive)
    [latsL, lonsL, namesL] = extract_latlon(lockers);
    [latsW, lonsW, namesW] = extract_latlon(warehouses);

    nL = numel(latsL);
    nW = numel(latsW);

    % --------- Pairwise distances (haversine) ----------
    % Output D(i,j) = distance from locker i to warehouse j (km)
    D = zeros(nL, nW);
    for i = 1:nL
        % vectorized across warehouses
        D(i, :) = haversine_km(latsL(i), lonsL(i), latsW(:), lonsW(:));
    end

    % Nearest warehouse for each locker
    [minDist, idxNearest] = min(D, [], 2);

    % Save distance matrix (row/col labels)
    outT = array2table(D, 'VariableNames', matlab.lang.makeValidName(namesW), ...
                       'RowNames', matlab.lang.makeValidName(namesL));
    try
        writetable(outT, 'distance_matrix.csv', 'WriteRowNames', true);
    catch
        warning('Could not write CSV to current folder.');
    end

    % --------- Plotting ----------
    make_plot(latsL, lonsL, namesL, latsW, lonsW, namesW, idxNearest, minDist);

    % --------- Return results ----------
    results.D = D;
    results.warehouses_coords = [latsW, lonsW];
    results.lockers_coords = [latsL, lonsL];
    results.lockers = lockers;
    results.warehouses = warehouses;
    results.namesLockers = namesL;
    results.namesWarehouses = namesW;
    results.nearestIdx = idxNearest;
    results.nearestDistKm = minDist;

    % also attempt to extract capacities and occupancy rates (if present)
    results.capacities = []; results.occupancyRate = [];
    try
        lowerVars = lower(lockers.Properties.VariableNames);
        capIdx = find(contains(lowerVars, {'capacity','cap'}), 1);
        occIdx = find(contains(lowerVars, {'occupancyrate','occupancy_rate','occupancyrate','occupancy','utilization','fillrate'}), 1);

        if ~isempty(capIdx)
            capRaw = lockers{:, capIdx};
            if iscell(capRaw), capRaw = str2double(capRaw); end
            results.capacities = double(capRaw);
        end
        if ~isempty(occIdx)
            occRaw = lockers{:, occIdx};
            if iscell(occRaw), occRaw = str2double(occRaw); end
            results.occupancyRate = double(occRaw);
        end
    catch
        % silently ignore extraction errors for these optional fields
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Helper functions
    function T = read_table_flexible(fname)
        % Read Excel or text files with common delimiters (comma, semicolon, tab).
        [~,~,ext] = fileparts(fname);
        if isempty(ext)
            ext = '.csv';
        end
        ext = lower(ext);
        try
            if any(strcmpi(ext, {'.xls','.xlsx','.xlsm','.xlsb'}))
                % Excel - let readtable handle it
                T = readtable(fname);
            else
                % For text files: try readtable first (auto-detect delimiter).
                try
                    T = readtable(fname);
                    return;
                catch
                    % if auto-detect fails, we'll try specific delimiters below
                end

                % Try a list of common delimiters
                delimiters = {',', ';', '\t'};
                success = false;
                lastErr = [];
                for k = 1:numel(delimiters)
                    try
                        opts = detectImportOptions(fname, 'FileType', 'text', 'Delimiter', delimiters{k});
                        % Some files may have no variable names in the first row; let readtable use opts
                        T = readtable(fname, opts);
                        success = true;
                        break;
                    catch ME
                        lastErr = ME;
                    end
                end

                if ~success
                    % As a last resort try reading as whitespace-delimited
                    try
                        opts = detectImportOptions(fname, 'FileType', 'text');
                        T = readtable(fname, opts);
                    catch ME
                        error('Failed to read text file "%s": %s', fname, ME.message);
                    end
                end
            end
        catch ME
            error('Failed to read file "%s": %s', fname, ME.message);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [lat, lon, ids] = extract_latlon(T)
        % find lat/long column names case-insensitively
        varNames = T.Properties.VariableNames;
        lowerVars = lower(varNames);
        % possible names
        latIdx = find(contains(lowerVars, 'lat'), 1);
        lonIdx = find(contains(lowerVars, {'lon','lng','long','longitude'}), 1);
        idIdx  = find(contains(lowerVars, {'id','name','locker_id','lockerid'}), 1);

        if isempty(latIdx) || isempty(lonIdx)
            error('Could not find lat/long columns in the file. Variable names found: %s', strjoin(varNames,', '));
        end

        lat = T{:, latIdx};
        lon = T{:, lonIdx};

        % If stored as strings, convert to numeric
        if iscell(lat), lat = str2double(lat); end
        if iscell(lon), lon = str2double(lon); end

        if any(isnan(lat)) || any(isnan(lon))
            warning('Some lat/lon values are NaN after conversion.');
        end

        if ~isempty(idIdx)
            idsRaw = T{:, idIdx};
            if isnumeric(idsRaw)
                ids = cellstr(num2str(idsRaw));
            else
                ids = cellstr(string(idsRaw));
            end
        else
            % create default names
            N = numel(lat);
            ids = arrayfun(@(k) sprintf('row_%d',k), (1:N)', 'UniformOutput', false);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function d = haversine_km(lat1, lon1, lat2, lon2)
        % Returns great-circle distance (km) between points
        % inputs can be scalars or vectors; output matches broadcasting similar to:
        % if lat1 scalar and lat2 vector -> returns vector
        R = 6371.0; % Earth's radius (km)

        % Ensure arrays are columns
        lat1 = double(lat1(:));
        lon1 = double(lon1(:));
        lat2 = double(lat2(:));
        lon2 = double(lon2(:));

        % If lat1 is scalar and lat2 is vector, broadcast lat1
        if isscalar(lat1) && ~isscalar(lat2)
            lat1 = repmat(lat1, size(lat2));
            lon1 = repmat(lon1, size(lon2));
        end
        if ~isscalar(lat1) && isscalar(lat2)
            lat2 = repmat(lat2, size(lat1));
            lon2 = repmat(lon2, size(lon1));
        end

        % If vectors of different orientation, make them match elementwise
        % Most callers will pass scalar/ vector combos or vectors of same length.
        dLat = deg2rad(lat2 - lat1);
        dLon = deg2rad(lon2 - lon1);
        a = sin(dLat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dLon/2).^2;
        c = 2 .* atan2(sqrt(a), sqrt(1-a));
        d = R .* c;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function make_plot(latsL, lonsL, namesL, latsW, lonsW, namesW, idxNearest, minDist)
        figure('Name','Lockers and Warehouses Map','NumberTitle','off');
        % If geoscatter + geobasemap available (Mapping Toolbox, R2019b+), use it
        if exist('geoscatter','file') && exist('geobasemap','file')
            geoscatter(latsL, lonsL, 5, 'filled', 'MarkerFaceColor',[0.85 0.33 0.10]); hold on;
            geoscatter(latsW, lonsW, 80, '^', 'filled', 'MarkerFaceColor',[0 0.45 0.74]);
            geobasemap('streets'); % or 'topographic','satellite' etc.
            legend('Lockers','Warehouses','Location','best');
            % label some points
            for m=1:numel(namesL)
                text(lonsL(m), latsL(m), namesL{m}, 'FontSize',8, 'HorizontalAlignment','left');
            end
            title('Lockers (orange) and Warehouses (blue)');
        else
            % Fallback: simple scatter with lat vs lon
            scatter(lonsL, latsL, 5, 'o','MarkerFaceColor',[0.85 0.33 0.10]); hold on;
            scatter(lonsW, latsW, 80, '^','MarkerFaceColor',[0 0.45 0.74]);
            for j=1:numel(latsL)
                plot([lonsL(j) lonsW(idxNearest(j))],[latsL(j) latsW(idxNearest(j))],'k-');
            end
            xlabel('Longitude'); ylabel('Latitude'); axis equal;
            legend('Lockers','Warehouses','Location','best');
            xlim([min([lonsL; lonsW]) - 0.01, max([lonsL; lonsW]) + 0.01]);
            ylim([min([latsL; latsW]) - 0.01, max([latsL; latsW]) + 0.01]);
            % add text labels
            for k=1:numel(namesL)
                text(lonsL(k)+0.0005, latsL(k)+0.0005, namesL{k}, 'FontSize',8);
            end
            for l=1:numel(namesW)
                text(lonsW(l)+0.0005, latsW(l)+0.0005, namesW{l}, 'FontSize',8, 'FontWeight','bold');
            end
            title('Lockers and Warehouses (lat/long scatter)');
            grid on;
            % Note: you can overlay an OSM tile image as background if you download one.
        end
    end

end