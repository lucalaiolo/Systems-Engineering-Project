classdef WarehouseServiceCompletionEvent < event
    %WAREHOUSESERVICECOMPLETIONEVENT Completion of a warehouse service activity.
    %   This event is scheduled when the picking/packing/loading process for an
    %   order finishes at a warehouse. The completion makes the order available
    %   for loading on a vehicle and potentially frees the warehouse server so
    %   that the next waiting order can begin service.

    properties
        warehouseId   % identifier of the warehouse where service completes
        orderIndex    % index of the order that just completed service
    end

    methods
        function obj = WarehouseServiceCompletionEvent(warehouseId, orderIndex, eventTime)
            if nargin < 3
                error('WarehouseServiceCompletionEvent requires warehouseId, orderIndex and time.');
            end
            label = sprintf('WarehouseServiceCompletion_W%d_O%d', warehouseId, orderIndex);
            parameters.warehouseId = warehouseId;
            parameters.orderIndex = orderIndex;
            obj@event(label, parameters, eventTime);
            obj.warehouseId = warehouseId;
            obj.orderIndex = orderIndex;
        end

        function randomVar = rnd(~)
            % The completion time is determined when the event is scheduled
            % (service time sampled at start of service). No additional random
            % variable is generated here.
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            simEngine.handleWarehouseServiceCompletion(obj.warehouseId, obj.orderIndex);
        end
    end
end