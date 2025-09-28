classdef futureEventsList < handle
    
    properties (Access = private, Constant)
        DEFAULT_CAPACITY = 16;
    end

    properties (Access = private)
        events_heap      % Cell array of event handles
        eventClocks      % Cached clock values for the heap
        heapsize         % Number of valid elements in the heap
        capacity         % Allocated capacity for the heap arrays
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = futureEventsList()
            obj.capacity = obj.DEFAULT_CAPACITY;
            obj.events_heap = cell(obj.capacity, 1);
            obj.eventClocks = zeros(obj.capacity, 1);
            obj.heapsize = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function that tells whether the heap is empty or not
        function isempty = Empty(obj)
            isempty = (obj.heapsize == 0);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Return the next event (without removing it from the heap)
        function nextEvent = First(obj)
            if obj.Empty()
                error('Future events list empty. There is no next event.');
            end
            nextEvent = obj.events_heap{1};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add a new event to the heap
        function Enqueue(obj, newEvent)
            if ~isa(newEvent, 'event')
                error('Future events list can only store objects derived from event.');
            end

            if obj.heapsize == obj.capacity
                obj.grow();
            end

            obj.heapsize = obj.heapsize + 1;
            obj.events_heap{obj.heapsize} = newEvent;
            obj.eventClocks(obj.heapsize) = newEvent.clock;
            obj.siftUp(obj.heapsize);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Dequeue an event
        function nextEvent = Dequeue(obj)
            if obj.Empty()
                error('Future events list empty. There is no next event.');
            end

            nextEvent = obj.events_heap{1};

            if obj.heapsize == 1
                obj.events_heap{1} = [];
                obj.eventClocks(1) = 0;
                obj.heapsize = 0;
                return;
            end

            obj.events_heap{1} = obj.events_heap{obj.heapsize};
            obj.eventClocks(1) = obj.eventClocks(obj.heapsize);
            obj.events_heap{obj.heapsize} = [];
            obj.eventClocks(obj.heapsize) = 0;
            obj.heapsize = obj.heapsize - 1;
            obj.siftDown(1);

            if obj.capacity > obj.DEFAULT_CAPACITY && obj.heapsize <= obj.capacity / 4
                obj.shrink();
            end
        end

    end %end methods

    methods (Access = private)
        
        function siftUp(obj, index)
            while index > 1
                parent = floor(index / 2);
                if obj.eventClocks(index) >= obj.eventClocks(parent)
                    break;
                end
                obj.swap(index, parent);
                index = parent;
            end
        end %end siftUp

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function siftDown(obj, index)
            while true
                left = 2 * index;
                right = left + 1;
                smallest = index;

                if left <= obj.heapsize && obj.eventClocks(left) < obj.eventClocks(smallest)
                    smallest = left;
                end

                if right <= obj.heapsize && obj.eventClocks(right) < obj.eventClocks(smallest)
                    smallest = right;
                end

                if smallest == index
                    break;
                end

                obj.swap(index, smallest);
                index = smallest;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function swap(obj, i, j)
            tmpEvent = obj.events_heap{i};
            obj.events_heap{i} = obj.events_heap{j};
            obj.events_heap{j} = tmpEvent;

            tmpClock = obj.eventClocks(i);
            obj.eventClocks(i) = obj.eventClocks(j);
            obj.eventClocks(j) = tmpClock;
        end %end swap

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function grow(obj)
            newCapacity = max(1, obj.capacity * 2);
            obj.events_heap = [obj.events_heap; cell(newCapacity - obj.capacity, 1)];
            obj.eventClocks = [obj.eventClocks; zeros(newCapacity - obj.capacity, 1)];
            obj.capacity = newCapacity;
        end %end grow

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function shrink(obj)
            proposed = max(obj.DEFAULT_CAPACITY, max(obj.heapsize, floor(obj.capacity / 2)));

            if proposed >= obj.capacity
                return;
            end

            obj.events_heap = obj.events_heap(1:proposed);
            obj.eventClocks = obj.eventClocks(1:proposed);
            obj.capacity = proposed;
        end %end shrink
    end %end private methods
end %end class definition

