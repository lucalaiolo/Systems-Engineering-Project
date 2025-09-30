classdef LockerDeliveryEvent < event
    %LOCKERDELIVERYEVENT Delivery of a batch of parcels to a locker stop.
    %   Each vehicle trip produces one LockerDeliveryEvent per visited locker.
    %   The event updates locker occupancy, records lead times and collects
    %   blocking statistics whenever the locker is full and parcels must be
    %   returned to the warehouse.

    properties
        warehouseId
        vehicleId
        stopIndex
    end

    methods
        function obj = LockerDeliveryEvent(warehouseId, vehicleId,stopIndex, eventTime)
            if nargin < 4
                error('LockerDeliveryEvent requires warehouseId, vehicleId, stopIndex and time.');
            end
            label = sprintf('LockerDelivery_W%d_V%d_Stop%d', warehouseId, vehicleId, stopIndex);
            params.warehouseId = warehouseId;
            params.vehicleId = vehicleId;
            params.stopIndex = stopIndex;
            obj@event(label, params, eventTime);
            obj.warehouseId = warehouseId;
            obj.vehicleId = vehicleId;
            obj.stopIndex = stopIndex;
        end

        function randomVar = rnd(~)
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            simEngine.handleLockerDelivery(obj.warehouseId, obj.vehicleId, obj.stopIndex);
        end
    end
end

