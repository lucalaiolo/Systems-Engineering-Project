classdef state < dynamicprops
    % state class
    % This class will represent the state of a system in the context of a
    % discrete event simulation.
    % It is a subclass of the class dynamicprops, which allows the properties 
    % to be defined and added during the call to the constructor

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = state(DataStruct)
            % Check whether input argument is valid or not
            if ~isstruct(DataStruct)
                error(['Error during the construction of a state object. ' ...
                    'The argument provided must be a struct.'])
            end
            fields = fieldnames(DataStruct);
            for i=1:length(fields)
                addprop(obj, fields{i});
                obj.(fields{i}) = DataStruct.(fields{i});
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update the state
        % This method can be used to clear the state properties too
        function update(obj, fields, newValues)
            % Check whether the inputs are valid or not
            % fields and newValues must be 1xn cell arrays
            if ~iscell(newValues) || ~isrow(newValues)
                error('Error during the update of the state. The newValues argument is not valid.')
            end
            if ~iscell(fields) || ~isrow(fields)
                error('Error during the update of the state. The fields argument is not valid.')
            end
            if length(fields) ~= length(newValues)
                error(['Error during the update of the state. The fields and newValues ', ...
                    'arguments do not have the same length']);
            end
            % Do the update
            for i=1:length(fields)
                field = fields{i};
                if ~isprop(obj, field)
                    error(['Error during the update of the state. One of the ' ...
                        'properties provided was not recognized.'])
                end
                obj.(field) = newValues{i};
            end
        end %end update

    end %end methods

end %end classdef

