classdef VehicleReturnEvent < event
    %VEHICLERETURNEVENT Event triggered when a vehicle completes its route and
    %returns to the warehouse. The event finalises the trip, reinserts blocked
    %parcels into the warehouse buffer (if any) and frees the vehicle so future
    %dispatches can occur.

    properties
        warehouseId
        vehicleId
    end

    methods
        function obj = VehicleReturnEvent(warehouseId, vehicleId, eventTime)
            if nargin < 3
                error('VehicleReturnEvent requires warehouseId, vehicleId and event time.');
            end
            label = sprintf('VehicleReturn_W%d_V%d', warehouseId, vehicleId);
            params.warehouseId = warehouseId;
            params.vehicleId = vehicleId;
            obj@event(label, params, eventTime);
            obj.warehouseId = warehouseId;
            obj.vehicleId = vehicleId;
        end

        function randomVar = rnd(~)
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            simEngine.handleVehicleReturn(obj.warehouseId, obj.vehicleId);
        end
    end
end
