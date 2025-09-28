function [orders, lockers, T_LW_minutes, T_LL_minutes] = read_and_preprocess(lockersFile, warehouseCoordsFile, ...
    ordersFile, speed)

lockers = readtable("parcel_lockers_stockholm_2.csv");

data = locker_warehouse_distances(lockersFile, warehouseCoordsFile);
numWarehouses = size(data.warehouses_coords);
numWarehouses = numWarehouses(1);
[lockerClustersIdx, C] = clusterizeLockers(data.lockers_coords, numWarehouses);
[~, locker2WarehouseAssignment] = assignClustersToWarehouses(C, data.warehouses_coords, ...
    lockerClustersIdx); 
data.lockers{:, "Warehouse_assigned"} = locker2WarehouseAssignment;
orders = read_orders_to_table_first(ordersFile);
orders = preprocessOrders(orders, locker2WarehouseAssignment);

[~, T_LW_minutes] = locker_warehouse_time('distance_matrix.csv', speed);
T_LW_minutes = table2array(T_LW_minutes);

[~, ~, T_LL_minutes] = locker_locker_time(data, speed);

end

