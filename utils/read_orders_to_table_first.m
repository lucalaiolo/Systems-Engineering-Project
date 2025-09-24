function [T, sortIdx] = read_orders_to_table_first(fname)
%READ_ORDERS_TO_TABLE_FIRST Read orders file and return table sorted by minutes since first timestamp.
%
%   [T, sortIdx] = read_orders_to_table_first(fname)
%
%   - Reads CSV/XLSX (tries readtable then common delimiters).
%   - Finds a timestamp-like column (case-insensitive).
%   - Parses timestamp strings into datetime (tries common formats then infers).
%   - Computes TimestampMinutes = minutes since earliest timestamp in file.
%   - Appends TimestampDatetime and TimestampMinutes to the table.
%   - Sorts rows by TimestampMinutes (ascending).
%   - Returns sorted table T and optional sortIdx such that T = T_original(sortIdx, :).

    arguments
        fname (1,:) char
    end

    if ~isfile(fname)
        error('File not found: %s', fname);
    end

    % --- read table (try readtable then fallbacks) ---
    try
        T = readtable(fname);
    catch
        % try common delimiters for text files
        delims = {',',';','\t'};
        readOK = false;
        lastErr = [];
        for k = 1:numel(delims)
            try
                opts = detectImportOptions(fname, 'FileType', 'text', 'Delimiter', delims{k});
                T = readtable(fname, opts);
                readOK = true;
                break;
            catch ME
                lastErr = ME;
            end
        end
        if ~readOK
            % last resort: try readtable without options (may error)
            try
                T = readtable(fname);
            catch ME
                error('Could not read file "%s": %s', fname, lastErr.message);
            end
        end
    end

    % --- find timestamp-like column (case-insensitive) ---
    vars = T.Properties.VariableNames;
    lowerVars = lower(vars);
    tIdx = find(contains(lowerVars, 'timestamp') | contains(lowerVars,'time') | strcmp(lowerVars,'ts'), 1);
    if isempty(tIdx)
        error('Could not find a Timestamp column. Available variables: %s', strjoin(vars,', '));
    end
    tsName = vars{tIdx};

    % --- get raw timestamp column and convert to string array ---
    rawTs = T{:, tIdx};
    if iscell(rawTs) || isstring(rawTs) || ischar(rawTs)
        rawStr = string(rawTs);
    elseif isnumeric(rawTs)
        % treat as datenum-like numeric
        try
            dt_numeric = datetime(rawTs, 'ConvertFrom','datenum');
            rawStr = string(dt_numeric);
        catch
            error('Unrecognized numeric timestamp format in column "%s".', tsName);
        end
    else
        rawStr = string(rawTs);
    end

    % --- try parsing with a few common formats ---
    formats = { ...
        'yyyy-MM-dd HH:mm', ...
        'yyyy-MM-dd HH:mm:ss', ...
        'yyyy/MM/dd HH:mm', ...
        'yyyy/MM/dd HH:mm:ss', ...
        'yyyy-MM-dd''T''HH:mm', ...
        'yyyy-MM-dd''T''HH:mm:ss' ...
        };

    parsed = false;
    for k = 1:numel(formats)
        try
            dt = datetime(rawStr, 'InputFormat', formats{k});
            if ~all(isnat(dt))
                parsed = true;
                break;
            end
        catch
            % ignore and try next format
        end
    end

    if ~parsed
        % final fallback: let MATLAB infer format
        try
            dt = datetime(rawStr);
        catch
            error('Failed to parse timestamps in column "%s". Example value: %s', tsName, rawStr(1));
        end
    end

    % --- compute minutes since earliest timestamp (first method) ---
    validMask = ~isnat(dt);
    if ~any(validMask)
        error('No valid timestamps parsed in column "%s".', tsName);
    end
    t0 = min(dt(validMask));
    minutesSinceFirst = minutes(dt - t0);  % numeric (may be fractional if seconds present)

    % --- sort by minutesSinceFirst (ascending) and build sorted table ---
    [~, sortIdx] = sort(minutesSinceFirst);     % indices into original table
    T_sorted = T(sortIdx, :);

    % --- attach parsed datetime and minutes to sorted table (in sorted order) ---
    % Ensure we create unique variable names that don't clash with existing ones
    origVars = T_sorted.Properties.VariableNames;
    dtVar = matlab.lang.makeUniqueStrings('TimestampDatetime', origVars);
    % Update origVars (include dtVar) before making minutes var unique
    origVars = [origVars, dtVar];
    minVar = matlab.lang.makeUniqueStrings('TimestampMinutes', origVars);

    T_sorted.(dtVar) = dt(sortIdx);
    T_sorted.(minVar) = minutesSinceFirst(sortIdx);

    % Return
    T = T_sorted;

    numOrders = size(T);
    numOrders = numOrders(1);
    eps = (1:numOrders)'*1e-6;

    % ensure the last column's contents are numeric (convert if needed)
    col = T{:, end};       % this extracts the contents (not a table)
    if ~isnumeric(col)
        col = str2double(col);   % convert cellstr/string to numeric (NaN on failure)
    end
    
    T{:, end} = col + eps;  % write back numeric values

end