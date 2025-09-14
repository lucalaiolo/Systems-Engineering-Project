classdef statsManager < handle
    % statsManager class
    
    properties
        statistics          % struct that maps field names to Statistic objects
        clock               % simulation clock
    end
    
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = statsManager(statStruct)
            if ~isstruct(statStruct)
                error(['Error during construction of statsManager object. ' ...
                    'The argument provided is not a struct.'])
            end
            obj.statistics = statStruct;
            obj.clock = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update method
        function update(obj, field, value)
            if nargin < 3
                value = 0;
            end
            if ~isfield(obj.statistics, field)
                error(['Error during the update of statsManager object. ' ...
                    'The field argument was not recognized.'])
            end
            obj.statistics.(field).update(value, obj.clock);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Clear method
        function clear(obj)
            fields = fieldnames(obj.statistics);
            for i=1:length(fields)
                obj.statistics.(fields{i}).clear();
            end
            obj.clock = 0;
        end

    end %end methods
end %end classdef

