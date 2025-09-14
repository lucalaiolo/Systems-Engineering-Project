classdef (Abstract) Statistic < handle
    % Abstract Statistic class
    methods (Abstract)
        update(obj, value, clock)
        clear(obj)
        result = getResult(obj);
    end
    
end
