classdef SumStatistic < Statistic
    % Subclass of Statistic class. This class collect sum statistics.
    properties
        sum = 0
    end
    
    methods

        function update(obj, toAdd, ~)
            obj.sum = obj.sum + toAdd;
        end

        function clear(obj)
            obj.sum = 0;
        end

        function result = getResult(obj)
            result = obj.sum;
        end

    end
end

