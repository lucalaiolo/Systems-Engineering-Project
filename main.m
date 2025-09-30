%% clear workspace
close all
clear
clc
rng default

%% Add dependencies
addpath(genpath('MATDES'));
addpath(genpath('Last Mile Delivery'));
addpath("utils\")

%% Import and preprocess data
lockersFile = 'parcel_lockers_stockholm_2.csv';
warehouseCoordsFile = 'warehouses_locations.xlsx';
ordersFile = 'stockholm_orders_sept2025.csv';
SPEED = 85;

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

onTimePercentage = zeros(nConfigs, 1);

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
        'capacity', Inf, ...
        'dispatchInterval', 3 * 60);
    warehouses = warehouses(cfg,:);
    %% Assemble the modelConfig structure
    modelConfig.orders = orders;
    modelConfig.lockers = lockers;
    %modelConfig.lockers{:,4} = modelConfig.lockers{:, 4} * 10;
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
    sim = LastMileDelivery(modelConfig);
    stopping = struct('maxTime', Inf);
    sim.run(stopping);

    %% Post processing
    [nOrders, ~] = size(orders);
    inTime = 0;
    for j=1:nOrders
        if ((sim.orders(j).deliveryTime - sim.orders(j).arrivalTime) < ...
            deliveryWindow(sim.orders(j).warehouseId))
            inTime = inTime + 1;
        end
    end
    onTimePercentage(i) = inTime / nOrders;

end
