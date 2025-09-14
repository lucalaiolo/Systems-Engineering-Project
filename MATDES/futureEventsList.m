classdef futureEventsList < handle

    % Implements a priority queue using a binary heap to manage a list of 
    % event objects. Events are dequeued according to their simulation clocks,
    % so that the earliest event is always the first one processed.

    properties
        events_heap
        heapsize
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = futureEventsList()
            obj.events_heap = event.empty;
            obj.heapsize = 0;
        end %end constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function that tells whether the heap is empty or not
        function isempty = Empty(obj)
            if obj.heapsize == 0
                isempty = true;
            else
                isempty = false;
            end
        end %end Empty
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Return the next event (without removing it from the heap)
        function nextEvent = First(obj)
            if ~obj.Empty()
                nextEvent = obj.events_heap(1);
            else
                error('Future events list empty. There is no next event.')
            end
        end %end First
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add a new event to the heap
        function Enqueue(obj, newEvent)
            obj.events_heap(end+1) = newEvent;
            obj.heapsize = obj.heapsize + 1;
            obj.rearrangeHeap(obj.heapsize)
        end %end Enqueue
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Dequeue an event
        function nextEvent = Dequeue(obj)
            if ~obj.Empty()
                nextEvent = obj.events_heap(1);
                obj.events_heap(1) = obj.events_heap(obj.heapsize);
                obj.events_heap = obj.events_heap(1:end-1);
                obj.heapsize = obj.heapsize - 1;
                obj.rearrangeHeap(1);
            else
               error('Future events list empty. There is no next event.');
            end
        end %end Dequeue
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rearrange the heap after dequeueing an event
        function rearrangeHeap(obj, i) %obj.heap is a heap except for the i-th position
            while (i > 1 && obj.events_heap(i).clock < obj.events_heap(Parent(i)).clock)
                obj.swap(i, Parent(i));
                i = Parent(i);
            end %end while
            while (Left(i) <= obj.heapsize && i ~= BestParentSons(obj, i))
                best = BestParentSons(obj, i);
                obj.swap(i, best);
                i = best;
            end %end while
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Nested auxiliary functions
            function parent = Parent(i)
                parent = floor((i-2)/2) + 1;
            end %end Parent
            
            function left = Left(i)
                left = 2*i;
            end %end Left

            function j = BestParentSons(obj, i)
                j = Left(i);
                k = j;
                if (k+1 <= obj.heapsize)
                    k = k + 1;
                end
                if obj.events_heap(k).clock < obj.events_heap(j).clock
                    j = k;
                end
                if obj.events_heap(i).clock <= obj.events_heap(j).clock
                    j = i;
                end
            end %end BestParentSons
            
        end %end rearrangeHeap
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Swap two elements of events_heap
        function swap(obj, i, j)
            tmp = obj.events_heap(i);
            obj.events_heap(i) = obj.events_heap(j);
            obj.events_heap(j) = tmp;
        end %end swap
    
    end %end methods
    
end %end class definition

