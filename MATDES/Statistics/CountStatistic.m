classdef CountStatistic < Statistic
    % Subclass of Statistic class. This class collect count statistics.
    properties
        count = 0
    end
    
    methods

        function update(obj, ~, ~)
            obj.count = obj.count + 1;
        end

        function clear(obj)
            obj.count = 0;
        end

        function result = getResult(obj)
            result = obj.count;
        end

    end
end

