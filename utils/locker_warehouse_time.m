function [T_hours, T_minutes] = locker_warehouse_time(distCsvFile, SPEED, varargin)
%LOCKER_WAREHOUSE_TIME  Convert distance matrix (km) to time matrices.
%   [T_hours, T_minutes] = locker_warehouse_time('distance_matrix.csv', 30);

% -------- Parse args --------
p = inputParser;
p.addParameter('OutDir','', @(s)ischar(s)||isstring(s));
p.addParameter('Prefix','time_matrix', @(s)ischar(s)||isstring(s));
p.addParameter('WriteFiles', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
OutDir = string(p.Results.OutDir);
Prefix = string(p.Results.Prefix);
WriteFiles = p.Results.WriteFiles;

if OutDir == ""
    [fdir,~,~] = fileparts(distCsvFile);
    if fdir == "", fdir = pwd; end
    OutDir = string(fdir);
end

if ~isscalar(SPEED) || ~isnumeric(SPEED) || ~isfinite(SPEED) || SPEED<=0
    error('SPEED must be a positive scalar in km/h.');
end

% -------- Read distance matrix (km) --------
% Read with row names; then coerce all columns to double safely.
D = readtable(distCsvFile, 'ReadRowNames', true);

% Coerce non-numeric columns to numeric (handles cases where numbers are stored as text)
for v = 1:width(D)
    if ~isnumeric(D.(v))
        D.(v) = str2double(string(D.(v)));
    end
end

% Verify numeric
if any(any(~isfinite(D{:,:})))
    warning('Some distance entries are non-finite after conversion.');
end

% -------- Compute time matrices --------
hours   = D{:,:} ./ SPEED;   % hours = km / (km/h)
minutes = hours * 60;

% Wrap back into tables with same labels
T_hours   = array2table(hours,   'RowNames', D.Properties.RowNames, ...
                                 'VariableNames', D.Properties.VariableNames);
T_minutes = array2table(minutes, 'RowNames', D.Properties.RowNames, ...
                                 'VariableNames', D.Properties.VariableNames);

% -------- Write outputs (optional) --------
if WriteFiles
    sp = regexprep(sprintf('%g',SPEED), '[^\w.]','_');
    f_hours   = fullfile(OutDir, sprintf('%s_hours_speed%s.csv',   Prefix, sp));
    f_minutes = fullfile(OutDir, sprintf('%s_minutes_speed%s.csv', Prefix, sp));
    writetable(T_hours,   f_hours,   'WriteRowNames', true);
    writetable(T_minutes, f_minutes, 'WriteRowNames', true);
    fprintf('Wrote:\n  %s\n  %s\n', f_hours, f_minutes);
end
end
