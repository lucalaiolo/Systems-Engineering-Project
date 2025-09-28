classdef LockerDeliveryEvent < event
    %LOCKERDELIVERYEVENT Delivery of a batch of parcels to a locker stop.
    %   Each vehicle trip produces one LockerDeliveryEvent per visited locker.
    %   The event updates locker occupancy, records lead times and collects
    %   blocking statistics whenever the locker is full and parcels must be
    %   returned to the warehouse.

    properties
        warehouseId
        stopIndex
    end

    methods
        function obj = LockerDeliveryEvent(warehouseId, stopIndex, eventTime)
            if nargin < 3
                error('LockerDeliveryEvent requires warehouseId, stopIndex and time.');
            end
            label = sprintf('LockerDelivery_W%d_Stop%d', warehouseId, stopIndex);
            params.warehouseId = warehouseId;
            params.stopIndex = stopIndex;
            obj@event(label, params, eventTime);
            obj.warehouseId = warehouseId;
            obj.stopIndex = stopIndex;
        end

        function randomVar = rnd(~)
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            simEngine.handleLockerDelivery(obj.warehouseId, obj.stopIndex);
        end
    end
end

