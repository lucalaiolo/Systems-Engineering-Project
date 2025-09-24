%%
clear
clc
close all
%%
addpath("MATDES\")
addpath("utils\")
%%
data = locker_warehouse_distances('parcel_lockers_stockholm_2.csv', 'base_config.xlsx');
%%
numWarehouses = size(data.warehouses_coords);
numWarehouses = numWarehouses(1);
[lockerClustersIdx, C] = clusterizeLockers(data.lockers_coords, numWarehouses);
[assignIdx, locker2WarehouseAssignment] = assignClustersToWarehouses(C, data.warehouses_coords, ...
    lockerClustersIdx); 
data.lockers{:, "Warehouse_assigned"} = locker2WarehouseAssignment;
%%
orders = read_orders_to_table_first("stockholm_orders_sept2025.csv");
orders = preprocessOrders(orders, locker2WarehouseAssignment);
%%
% Get time distance matrix
SPEED = 65; % speed of the vehicles, in kilometers per hour
[T_LW_hours, T_LW_minutes] = locker_warehouse_time('distance_matrix.csv', SPEED);

[D_LL, T_LL_hours, T_LL_minutes] = locker_locker_time(data, SPEED);
