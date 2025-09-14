classdef MinMaxStatistic < Statistic
    % Subclass of Statistic class. This class collect minimum/maximum statistics.
    properties
        minimum = Inf
        maximum = -Inf 
    end
    
    methods

        function update(obj, value, ~)
            obj.minimum = min(obj.minimum, value);
            obj.maximum = max(obj.maximum, value);
        end

        function clear(obj)
            obj.minimum = Inf;
            obj.maximum = -Inf;
        end

        function result = getResult(obj)
            result = [obj.minimum, obj.maximum];
        end

    end
end
