function [D_LL, T_LL_hours, T_LL_minutes] = locker_locker_time(data, SPEED, varargin)
%LOCKER_LOCKER_TIME Pairwise distances and times between lockers.
% Inputs:
%   data.lockers_coords : [n x 2] matrix [lat, lon] in degrees
%   data.namesLockers   : (optional) cellstr/string array of locker IDs
%   SPEED               : scalar km/h
% Name-Value:
%   OutDir   : output folder (default: pwd)
%   Prefix   : filename prefix (default: 'locker_lockers')
%   WriteCSV : true/false (default: true)
%
% Outputs:
%   D_LL         : [n x n] distance matrix (km)
%   T_LL_hours   : [n x n] time matrix (hours)
%   T_LL_minutes : [n x n] time matrix (minutes)

p = inputParser;
p.addParameter('OutDir', pwd, @(s)ischar(s)||isstring(s));
p.addParameter('Prefix','locker_lockers', @(s)ischar(s)||isstring(s));
p.addParameter('WriteCSV', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
outDir   = string(p.Results.OutDir);
prefix   = string(p.Results.Prefix);
writeCSV = p.Results.WriteCSV;

% --- Inputs & checks ---
coords = data.lockers_coords;
if size(coords,2) ~= 2
    error('data.lockers_coords must be n x 2 [lat, lon].');
end
lat = double(coords(:,1));
lon = double(coords(:,2));
n   = numel(lat);

if ~isscalar(SPEED) || ~isfinite(SPEED) || SPEED <= 0
    error('SPEED must be a positive scalar in km/h.');
end

% Locker names (optional)
if isfield(data,'namesLockers') && ~isempty(data.namesLockers)
    names = cellstr(string(data.namesLockers));
else
    names = arrayfun(@(k)sprintf('L%03d',k), (1:n)', 'UniformOutput', false);
end
% Valid, unique MATLAB variable names for columns
varNames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names));

% --- Pairwise haversine (vectorized) ---
R = 6371.0;                 % km
phi = deg2rad(lat(:));      % n x 1
lam = deg2rad(lon(:));

dphi = phi - phi.';         % n x n
dlam = lam - lam.';

a = sin(dphi/2).^2 + cos(phi) * cos(phi.').* sin(dlam/2).^2;
c = 2*atan2(sqrt(a), sqrt(1-a));
D_LL = R * c;
D_LL(1:n+1:end) = 0;        % exact zeros on diagonal

% --- Times ---
T_LL_hours   = D_LL ./ SPEED;
T_LL_minutes = T_LL_hours * 60;

% --- Save labeled CSVs (optional) ---
if writeCSV
    T_D  = array2table(D_LL,       'RowNames', names, 'VariableNames', varNames);
    T_H  = array2table(T_LL_hours, 'RowNames', names, 'VariableNames', varNames);
    T_M  = array2table(T_LL_minutes,'RowNames', names, 'VariableNames', varNames);

    fD = fullfile(outDir, sprintf('%s_distance_km.csv', prefix));
    fH = fullfile(outDir, sprintf('%s_time_hours_speed%g.csv', prefix, SPEED));
    fM = fullfile(outDir, sprintf('%s_time_minutes_speed%g.csv', prefix, SPEED));

    writetable(T_D, fD, 'WriteRowNames', true);
    writetable(T_H, fH, 'WriteRowNames', true);
    writetable(T_M, fM, 'WriteRowNames', true);
    fprintf('Wrote:\n  %s\n  %s\n  %s\n', fD, fH, fM);
end
end
