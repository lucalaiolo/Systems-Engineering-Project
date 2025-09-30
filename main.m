%% clear workspace
close all
clear
clc
rng default

%% Add dependencies
addpath(genpath('MATDES'));
addpath(genpath('Last Mile Delivery'));
addpath(genpath('utils'));

%% Import and preprocess data
lockersFile = 'parcel_lockers_stockholm_2.csv';
warehouseCoordsFile = 'warehouses_locations.xlsx';
ordersFile = 'stockholm_orders_sept2025.csv';
SPEED = 40;

warehousesCoordinates = readtable('warehouses_locations.xlsx');
[numWarehouses, ~] = size(warehousesCoordinates);
configs = {};
configs{end+1} = [1, 2]';
for i=3:numWarehouses
    configs{end+1} = [1, 2, i]';
end
for i=3:numWarehouses
    for j=(i+1):numWarehouses
        configs{end+1} = [1,2,i,j]';
    end
end

nConfigs = length(configs);
existingWData = table2array(readtable("warehouse_data__Existing_Warehouses.csv"));
newWData = table2array(readtable("warehouse_data__New_Warehouses.csv"));

numExistingWarehouses = size(existingWData, 1);

configCosts = zeros(nConfigs, 1);
for idx = 1:nConfigs
    configCosts(idx) = computeConfigurationCost(configs{idx}, existingWData, newWData);
end

nReplications = 2; % Number of simulation replications per configuration
onTimePercentage = nan(nConfigs, nReplications);
simulationResults = repmat(struct( ...
    'config', [], ...
    'onTimePercentage', nan(1, nReplications), ...
    'meanOnTimePercentage', NaN, ...
    'stdOnTimePercentage', NaN, ...
    'meanLockerBlockingProbability', [], ...
    'cost', NaN, ...
    'lockers', []), nConfigs, 1);

lockerBlockingProbabilities = cell(nConfigs, 1);

for idx=1:nConfigs
    simulationResults(idx).config = configs{idx};
    simulationResults(idx).cost = configCosts(idx);
end

for i=1:nConfigs
    cfg = configs{i};
    deliveryWindow = zeros(length(cfg), 1); 
    for j=1:length(cfg)
        if cfg(j) < 3
            deliveryWindow(j) = existingWData(cfg(j), 4);
        else
            deliveryWindow(j) = newWData(cfg(j)-2, 7);
        end
    end
    
    deliveryWindow = deliveryWindow * 60;

    [orders, lockers, T_LW_minutes, T_LL_minutes] = read_and_preprocess(lockersFile, ...
        warehouseCoordsFile, ordersFile, SPEED, cfg);

    warehouses = readtable('warehouses_locations.xlsx');
    warehouses{:, 'vehicle'} = struct( ...
        'capacity', 100, ...
        'dispatchInterval', 60, ...
        'fleetSize', 10);
    warehouses = warehouses(cfg,:);
    %% Assemble the modelConfig structure
    modelConfig.orders = orders;
    modelConfig.lockers = lockers;
    modelConfig.warehouses = warehouses;
    
    travelTimes.warehouseToLocker = T_LW_minutes';
    travelTimes.lockerToLocker = T_LL_minutes;
    travelTimes.lockerToWarehouse = T_LW_minutes;
    
    modelConfig.travelTimes = travelTimes;
    
    modelConfig.serviceTimesFiles = {'warehouse_data__Existing_Warehouses.csv',...
        'warehouse_data__New_Warehouses.csv'};
    
    modelConfig.wCfg = cfg;

    modelConfig.options = struct( ...
        'storeLeadTimes', true, ...
        'stopWhenAllDelivered', true);
    
    [nOrders, ~] = size(orders);
    numLockers = height(lockers);
    if isempty(lockerBlockingProbabilities{i})
        lockerBlockingProbabilities{i} = nan(nReplications, numLockers);
    end

    for replication=1:nReplications
        sim = LastMileDelivery(modelConfig);
        stopping = struct('maxTime', Inf);
        sim.run(stopping);

        lockerStats = sim.getLockerBlockingSummary();
        if isfield(lockerStats, 'probability') && ~isempty(lockerStats.probability)
            lockerBlockingProbabilities{i}(replication, :) = lockerStats.probability(:)';
        end

        inTime = 0;
        for j=1:nOrders
            if ((sim.orders(j).deliveryTime - sim.orders(j).arrivalTime) < ...
                deliveryWindow(sim.orders(j).warehouseId))
                inTime = inTime + 1;
            end
        end
        onTimePercentage(i, replication) = inTime / nOrders;
    end

    simulationResults(i).config = cfg;
    simulationResults(i).onTimePercentage = onTimePercentage(i, :);
    simulationResults(i).meanOnTimePercentage = mean(onTimePercentage(i, :), 'omitnan');
    simulationResults(i).stdOnTimePercentage = std(onTimePercentage(i, :), 'omitnan');
    if ~isempty(lockerBlockingProbabilities{i})
        simulationResults(i).meanLockerBlockingProbability = mean(lockerBlockingProbabilities{i}, 1, 'omitnan');
    end
    simulationResults(i).lockers = lockers;
    
end

meanOnTimePercentage = arrayfun(@(result) result.meanOnTimePercentage, simulationResults);
stdOnTimePercentage = arrayfun(@(result) result.stdOnTimePercentage, simulationResults);
costPerConfiguration = arrayfun(@(result) result.cost, simulationResults);

configLabels = cellfun(@(cfg) mat2str(cfg'), configs, 'UniformOutput', false);
summaryTable = table(configLabels(:), meanOnTimePercentage(:), stdOnTimePercentage(:), ...
    costPerConfiguration(:), ...
    'VariableNames', {'Configuration', 'MeanOnTimePercentage', 'StdOnTimePercentage', 'Cost'});

validConfigurations = find(~isnan(meanOnTimePercentage));
if ~isempty(validConfigurations)
    [~, bestConfigLocalIdx] = max(meanOnTimePercentage(validConfigurations));
    bestConfigIdx = validConfigurations(bestConfigLocalIdx);
    bestBlockingProbabilities = simulationResults(bestConfigIdx).meanLockerBlockingProbability;
    lockersForPlot = simulationResults(bestConfigIdx).lockers;
    plotLockerBlockingMap(lockersForPlot, bestBlockingProbabilities);

    validCosts = costPerConfiguration(validConfigurations);
    validMeanOnTime = meanOnTimePercentage(validConfigurations);
    [sortedCosts, order] = sort(validCosts);
    sortedMeanOnTime = validMeanOnTime(order);
    bestOnTimeEnvelope = cummax(sortedMeanOnTime);

    figure;
    plot(sortedCosts, bestOnTimeEnvelope, '-o', 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'auto');
    grid on;
    xlabel('Total Cost (Capex + Opex)');
    ylabel('Best On-Time Delivery Percentage');
    title('Best Achievable On-Time Delivery by Total Cost');
else
    warning('No valid simulation results available to plot locker blocking map or cost-performance curve.');
end

%% Helper functions
function plotLockerBlockingMap(lockerTable, blockingProbabilities)
    if isempty(lockerTable) || isempty(blockingProbabilities)
        warning('Locker data or blocking probabilities are not available for plotting.');
        return;
    end
    
    validMask = ~isnan(blockingProbabilities);
    if ~any(validMask)
        warning('All locker blocking probabilities are NaN. Skipping map plot.');
        return;
    end
    
    scatterSizes = 60;
    figure;
    scatter(lockerTable.Longitude(validMask), lockerTable.Latitude(validMask), scatterSizes, ...
        blockingProbabilities(validMask), 'filled');
    colorbar;
    colormap(parula);
    maxBlocking = max(blockingProbabilities(validMask));
    if maxBlocking <= 0
        maxBlocking = 1;
    end
    caxis([0, maxBlocking]);
    title('Locker Blocking Probability Map');
    xlabel('Longitude');
    ylabel('Latitude');
    axis equal;
end