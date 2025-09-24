function orders = preprocessOrders(orders, lockerAssignment)

numStr = regexprep(orders.Locker_ID, '\D', '');   % e.g., 'PL008' -> '008'
lockerNum = str2double(numStr);
orders.Locker_ID = lockerNum;           % -> 8 (leading zeros handled)

% Valid indices into lockerAssignment
valid = ~isnan(lockerNum) & lockerNum >= 1 & lockerNum <= numel(lockerAssignment);

% Preallocate as double
orders.Warehouse = nan(height(orders), 1);

% Assign
orders.Warehouse(valid) = lockerAssignment(lockerNum(valid));

end

