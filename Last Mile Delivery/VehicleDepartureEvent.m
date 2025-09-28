classdef VehicleDepartureEvent < event
    %VEHICLEDEPARTUREEVENT Departure of a delivery vehicle from a warehouse.
    %   Vehicles can be dispatched because the load reached capacity or because
    %   a scheduled departure time elapsed. The trigger is stored in the
    %   triggerType property so that the LastMileDelivery model can apply the
    %   appropriate logic (e.g., rescheduling future departures when the event
    %   corresponds to a timetable trigger).

    properties
        warehouseId
        triggerType % 'schedule' or 'capacity'
    end

    methods
        function obj = VehicleDepartureEvent(warehouseId, triggerType, eventTime)
            if nargin < 3
                error('VehicleDepartureEvent requires warehouseId, triggerType and time.');
            end
            if nargin < 2 || isempty(triggerType)
                triggerType = 'schedule';
            end
            label = sprintf('VehicleDeparture_W%d_%s', warehouseId, triggerType);
            params.warehouseId = warehouseId;
            params.triggerType = triggerType;
            obj@event(label, params, eventTime);
            obj.warehouseId = warehouseId;
            obj.triggerType = triggerType;
        end

        function randomVar = rnd(~)
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            simEngine.handleVehicleDeparture(obj.warehouseId, obj.triggerType);
        end
    end
end
