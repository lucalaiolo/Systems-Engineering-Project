classdef SimulationEngine < handle
    % SimulationEngine
    % The SimulationEngine class defines the framework for discrete event 
    % simulation

    properties
        state                   % a state object
        clock = 0               % simulation clock
        eventsList              % a futureEventsList object
        stats                   % a statsManager object
    end
    
    methods (Abstract)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Termination condition
        flag = terminationCondition(obj, stoppingParameters)

    end

    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = SimulationEngine(stateStruct, statStruct)
            % Initialize the state
            if ~isstruct(stateStruct) || ~isstruct(statStruct)
                error(['Error during construction of SimulationEngine object. ' ...
                    'One of the arguments is not a struct.'])
            end
            obj.state = state(stateStruct);
            
            % Initialize the clock
            obj.clock = 0;
            
            % Initialize the future events list
            obj.eventsList = futureEventsList();
            
            % Initialize the statistics
            obj.stats = statsManager(statStruct);
            
        end %end constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialize the events list
        % events must be an array of subclasses of event class
        % This function must be called after the creation of a
        % SimulationEngine object and before calling the run method
        function initEventsList(obj, events)
            for i=1:length(events)
                obj.eventsList.Enqueue(events(i));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run the simulation
        function run(obj, stoppingParameters)
            while ~obj.terminationCondition(stoppingParameters)
                % Determine the next event
                nextEvent = obj.eventsList.Dequeue();
                % Advance the simulation clock
                obj.clock = nextEvent.clock;
                obj.stats.clock = obj.clock;
                % Update the system state and schedule the next event
                % The statistics will be update inside this function
                nextEvent.manageEvent(obj);
            end
            % Print a report
            obj.generateReport();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Clear the state and the statistical counters
        % The state will be cleared using the struct DataStruct
        function clear(obj, DataStruct)
            obj.clock = 0;
            obj.eventsList = futureEventsList();
            fields = fieldnames(DataStruct)';
            values = struct2cell(DataStruct)';
            obj.state.update(fields, values);
            obj.stats.clear();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate final report
        function generateReport(obj)
            disp('Simulation ended.')
            disp('Statistics collected:');
            tracked = fieldnames(obj.stats.statistics);
            for i=1:length(tracked)
                field = tracked{i};
                disp([field, ': ']);
                disp(num2str(obj.stats.statistics.(field).getResult()));
            end
        end

    end %end methods

end %end classdef

