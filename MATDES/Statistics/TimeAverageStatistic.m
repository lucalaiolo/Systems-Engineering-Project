classdef TimeAverageStatistic < Statistic
    % Subclass of Statistic class. This class collect time average statistics.
    properties
        weightedSum = 0
        totalTime = 0
        lastUpdateTime = 0
    end
    
    methods

        function update(obj, oldValue, clock)
            deltaT = clock - obj.lastUpdateTime;
            obj.weightedSum = obj.weightedSum + deltaT * oldValue;
            obj.totalTime = obj.totalTime + deltaT;
            obj.lastUpdateTime = clock;
        end

        function clear(obj)
            obj.weightedSum = 0;
            obj.totalTime = 0;
            obj.lastUpdateTime = 0;
        end

        function result = getResult(obj)
            if obj.totalTime == 0
                error('Division by zero encountered.')
            end
            result = obj.weightedSum / obj.totalTime;
        end

    end
end


