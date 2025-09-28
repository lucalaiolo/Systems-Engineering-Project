classdef VehicleReturnEvent < event
    %VEHICLERETURNEVENT Event triggered when a vehicle completes its route and
    %returns to the warehouse. The event finalises the trip, reinserts blocked
    %parcels into the warehouse buffer (if any) and frees the vehicle so future
    %dispatches can occur.

    properties
        warehouseId
    end

    methods
        function obj = VehicleReturnEvent(warehouseId, eventTime)
            if nargin < 2
                error('VehicleReturnEvent requires warehouseId and event time.');
            end
            label = sprintf('VehicleReturn_W%d', warehouseId);
            params.warehouseId = warehouseId;
            obj@event(label, params, eventTime);
            obj.warehouseId = warehouseId;
        end

        function randomVar = rnd(~)
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            simEngine.handleVehicleReturn(obj.warehouseId);
        end
    end
end
