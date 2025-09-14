classdef (Abstract) event < handle & matlab.mixin.Heterogeneous
    % event
    % event is a subclass of the matlab.mixin.Heterogeneous class because
    % this class supports forming heterogeneous arrays. A heterogeneous 
    % array is an array of objects that differ in their specific class, but
    % are all derived from or are instances of a root class.
    % In the context of discrete event simulation, it is useful to build
    % heterogeneous arrays of specific events classes derived from the event 
    % superclass (see the futureEventsList class for example).

    properties
        label               % label for the event (typically, a string)
        parameters          % parameters used to generate random variables
        clock               % event clock
    end
    
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = event(label, parameters, clock)
            if nargin < 3
                clock = 0;
            end
            obj.label = label;
            obj.parameters = parameters;
            obj.clock = clock;
        end
        
    end

    methods (Abstract)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function to generate random variables
        randomVar = rnd(obj)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function to manage the event updating a SimulationEngine object and
        % to schedule the next events
        manageEvent(obj, simEngine)

    end

end

