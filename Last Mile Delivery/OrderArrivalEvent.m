classdef OrderArrivalEvent < event
    %ORDERARRIVALEVENT Event representing the arrival of an order at the system.
    %   Orders are assumed to come from the preprocessed empirical stream and
    %   are scheduled deterministically at the timestamps contained in the
    %   configuration. Each event is responsible for inserting a single order
    %   into the warehouse queue and scheduling the next arrival if any order
    %   remains in the stream.

    properties
        orderIndex   % index of the order within the simulation data structure
    end

    methods
        function obj = OrderArrivalEvent(orderIndex, eventTime)
            % Constructor for the order arrival event.
            if nargin < 2
                error('OrderArrivalEvent requires the order index and event time.');
            end
            label = sprintf('OrderArrival_%d', orderIndex);
            obj@event(label, struct(), eventTime);
            obj.orderIndex = orderIndex;
        end

        function randomVar = rnd(~)
            % No random variable is generated directly by this event because
            % inter-arrival times are deterministic (taken from the data
            % stream). The SimulationEngine framework, however, expects a
            % method to be defined, so we return an empty array.
            randomVar = [];
        end

        function manageEvent(obj, simEngine)
            % Dispatch management of the event to the simulation model.
            simEngine.handleOrderArrival(obj.orderIndex);
        end
    end
end
