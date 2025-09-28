classdef LockerPickupEvent < event
    %LOCKERPICKUPEVENT Event representing a customer pickup completion.
    %   The event can optionally carry the order index associated with the
    %   pickup so downstream logic can update the specific parcel record. If
    %   the order index is omitted (legacy behaviour), the simulator will fall
    %   back to removing the first parcel present in the locker.

    properties
        lockerId
        orderIndex
    end

    methods
        function obj = LockerPickupEvent(lockerId, varargin)
            % Constructor supporting both the legacy signature
            %   LockerPickupEvent(lockerId, eventTime)
            % and the per-order signature
            %   LockerPickupEvent(lockerId, orderIndex, eventTime).

            if nargin < 2
                error('LockerPickupEvent requires at least lockerId and event time.');
            end

            if isscalar(varargin)
                % Legacy form: lockerId, eventTime
                orderIndex = 0;
                eventTime = varargin{1};
            elseif numel(varargin) >= 2
                orderIndex = varargin{1};
                eventTime = varargin{2};
            else
                error('LockerPickupEvent received an unsupported argument configuration.');
            end

            if ~isscalar(orderIndex) || ~isnumeric(orderIndex)
                error('orderIndex must be a numeric scalar.');
            end

            if orderIndex > 0
                label = sprintf('LockerPickup_L%d_O%d', lockerId, orderIndex);
            else
                label = sprintf('LockerPickup_L%d', lockerId);
            end

            params.lockerId = lockerId;
            params.orderIndex = orderIndex;
            obj@event(label, params, eventTime);
            obj.lockerId = lockerId;
            obj.orderIndex = orderIndex;
        end

        function randomVar = rnd(~)
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            simEngine.handleLockerPickup(obj.lockerId, obj.orderIndex);
        end
    end
end


