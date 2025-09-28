%% clear workspace
close all
clear
clc

%% Add dependencies
addpath(genpath('MATDES'));
addpath(genpath('Last Mile Delivery'));
addpath("utils\")

%% Import and preprocess data
lockersFile = 'parcel_lockers_stockholm_2.csv';
warehouseCoordsFile = 'base_config.xlsx';
ordersFile = 'stockholm_orders_sept2025.csv';
SPEED = 55;

[orders, lockers, T_LW_minutes, T_LL_minutes] = read_and_preprocess(lockersFile, ...
    warehouseCoordsFile, ordersFile, SPEED);

warehouses = readtable('base_config.xlsx');
warehouses{:, 'vehicle'} = struct( ...
    'capacity', Inf, ...
    'dispatchInterval', 10 * 60);

%% Assemble the modelConfig structure
modelConfig.orders = orders;
modelConfig.lockers = lockers;
modelConfig.warehouses = warehouses;

travelTimes.warehouseToLocker = T_LW_minutes';
travelTimes.lockerToLocker = T_LL_minutes;
travelTimes.lockerToWarehouse = T_LW_minutes;

modelConfig.travelTimes = travelTimes;

modelConfig.options = struct( ...
    'storeLeadTimes', true, ...
    'stopWhenAllDelivered', true);

%% Instantiate the simulation
sim = LastMileDelivery(modelConfig);
stopping = struct('maxTime', Inf);
rng default
sim.run(stopping);

%% Post processing
