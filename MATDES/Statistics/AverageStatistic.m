classdef AverageStatistic < Statistic
    % Subclass of Statistic class. This class collect average statistics.
    properties
        count = 0
        sum = 0
    end
    
    methods

        function update(obj, toAdd, ~)
            obj.count = obj.count + 1;
            obj.sum = obj.sum + toAdd;
        end

        function clear(obj)
            obj.count = 0;
            obj.sum = 0;
        end

        function result = getResult(obj)
            if obj.count == 0
                error('Division by zero encountered.')
            end
            result = obj.sum / obj.count;
        end

    end
end
