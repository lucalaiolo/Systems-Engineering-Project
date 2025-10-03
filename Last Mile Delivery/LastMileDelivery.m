classdef LastMileDelivery < SimulationEngine
    %LASTMILEDELIVERY Discrete-event simulator for the last-mile delivery system.
    %   The model leverages the MATDES library to represent the interaction
    %   between warehouses (infinite-server queues fed by an empirical order
    %   stream), capacitated vehicles operating batched delivery routes and
    %   parcel lockers with stochastic customer pick-ups. The class orchestrates
    %   the creation of events, maintains the system state and collects key
    %   performance indicators such as queue lengths, utilizations, lead times
    %   and blocking probabilities.
    
    properties
        config              % configuration data for warehouses and lockers
        orders              % struct array containing the orders to simulate
        travelTimes         % travel-time matrices (warehouse-locker, locker-locker)
        options             % simulation options (max time, storage behaviour, ...)
        statNames           % cached statistic field names for quick lookup
        leadTimeSamples     % vector of collected lead-time samples
        leadTimeCount       % number of stored lead-time samples
        completedOrders     % number of orders successfully delivered
        totalOrders         % total number of orders in the arrival stream
        blockedParcels      % number of parcels blocked at lockers
        randomStream        % optional RandStream used for reproducible randomness
    end

    methods

        function obj = LastMileDelivery(modelConfig)
            % Constructor. The input must be a struct with fields describing the
            % orders, the physical network and the stochastic assumptions used by
            % the simulation. The constructor validates the configuration,
            % initialises the MATDES SimulationEngine superclass and sets up the
            % initial event list.

            parsed = LastMileDelivery.parseModelConfig(modelConfig);
            stateStruct = LastMileDelivery.createInitialState(parsed);
            statStruct = LastMileDelivery.createStatistics(parsed);

            obj@SimulationEngine(stateStruct, statStruct);

            obj.config = parsed.config;
            obj.orders = parsed.orders;
            obj.travelTimes = parsed.travelTimes;
            obj.options = parsed.options;
            obj.randomStream = parsed.randomStream;
            [obj.totalOrders, ~] = size(obj.orders);
            obj.completedOrders = 0;
            obj.blockedParcels = 0;
            obj.leadTimeCount = 0;
            if obj.options.storeLeadTimes
                obj.leadTimeSamples = nan(max(1, obj.totalOrders), 1);
            else
                obj.leadTimeSamples = [];
            end

            obj.prepareStatisticNames();

            % Schedule the initial set of events (first order arrival, vehicle
            % departure timers, locker pick-ups if initial inventory is present).
            initialEvents = obj.createInitialEvents();
            if ~isempty(initialEvents)
                obj.initEventsList(initialEvents);
            end
        end %end constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function flag = terminationCondition(obj, stoppingParameters)
            %TERMINATIONCONDITION Stop when there are no future events or when
            %the provided stopping criteria are met (max time or max number of
            %delivered orders).

            if obj.eventsList.Empty()
                flag = true;
                return;
            end

            if nargin < 2 || isempty(stoppingParameters)
                stoppingParameters = struct();
            end

            flag = false;

            if isstruct(stoppingParameters)
                if isfield(stoppingParameters, 'maxTime') && ~isempty(stoppingParameters.maxTime)
                    nextEvent = obj.eventsList.First();
                    if nextEvent.clock > stoppingParameters.maxTime
                        flag = true;
                        return;
                    end
                end
                if isfield(stoppingParameters, 'maxDeliveries') && ~isempty(stoppingParameters.maxDeliveries)
                    if obj.completedOrders >= stoppingParameters.maxDeliveries
                        flag = true;
                        return;
                    end
                end
            end

            if obj.options.stopWhenAllDelivered && obj.completedOrders >= obj.totalOrders
                flag = true;
                return;
            end

        end %end terminationCondition
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generateReport(obj)
            %GENERATEREPORT Custom report summarising the main performance
            %indicators along with the statistics managed by MATDES.
            
            fprintf('Simulation ended.');
            fprintf('Total orders delivered: %d / %d', obj.completedOrders, obj.totalOrders);
            fprintf('Total blocked parcels returned: %d', obj.blockedParcels);

            if obj.options.storeLeadTimes && obj.leadTimeCount > 0
                sample = obj.leadTimeSamples(1:obj.leadTimeCount);
                avgLead = mean(sample);
                sorted = sort(sample);
                idx95 = max(1, ceil(0.95 * numel(sorted)));
                perc95 = sorted(idx95);
                fprintf('Average lead time (minutes): %.3f', avgLead);
                fprintf('95th percentile lead time (minutes): %.3f', perc95);
            end

            disp('Statistics collected:');
            tracked = fieldnames(obj.stats.statistics);
            for i = 1:numel(tracked)
                field = tracked{i};
                try
                    value = obj.stats.statistics.(field).getResult();
                    disp([field, ': ']);
                    disp(num2str(value));
                catch ME
                    warning('Could not compute statistic %s: %s', field, ME.message);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lockerStats = getLockerBlockingSummary(obj)
            %GETLOCKERBLOCKINGSUMMARY Return aggregated locker blocking data.
            %   The returned struct contains the total number of blocked parcel
            %   attempts, the total number of delivery attempts and the
            %   resulting blocking probability for each locker.

            numLockers = obj.config.numLockers;
            blockedAttempts = zeros(numLockers, 1);
            totalAttempts = zeros(numLockers, 1);
            probability = nan(numLockers, 1);

            for l = 1:numLockers
                field = obj.statNames.lockerBlocking{l};
                statObj = obj.stats.statistics.(field);
                blockedAttempts(l) = statObj.sum;
                totalAttempts(l) = statObj.count;
                if totalAttempts(l) > 0
                    probability(l) = blockedAttempts(l) / totalAttempts(l);
                end
            end

            lockerStats = struct('blockedAttempts', blockedAttempts, ...
                'totalAttempts', totalAttempts, 'probability', probability);
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function handleOrderArrival(obj, orderIndex)
            %HANDLEORDERARRIVAL Process an order arrival event.

            if orderIndex > obj.totalOrders || orderIndex < 1
                return;
            end

            order = obj.orders(orderIndex,:);
            warehouseId = order.warehouseId;

            % Update queue length time-average before modifying the queue.
            queueField = obj.statNames.warehouseQueueTimeAvg{warehouseId};
            oldQueueLength = obj.state.warehouseQueueLengths(warehouseId);
            obj.stats.update(queueField, oldQueueLength);

            % In an M/M/inf setting every order can begin service immediately,
            % therefore no waiting queue is maintained.
            obj.startWarehouseService(warehouseId, orderIndex);

            % Schedule the next arrival from the empirical stream.
            if orderIndex < obj.totalOrders
                nextTime = obj.orders(orderIndex + 1).arrivalTime;
                obj.eventsList.Enqueue(OrderArrivalEvent(orderIndex + 1, nextTime));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function handleWarehouseServiceCompletion(obj, warehouseId, orderIndex)
            %HANDLEWAREHOUSESERVICECOMPLETION Manage completion of a warehouse
            %service activity.

            if warehouseId < 1 || warehouseId > obj.config.numWarehouses
                return;
            end

            obj.state.warehouseInService(warehouseId) = max(0, obj.state.warehouseInService(warehouseId) - 1);
            obj.orders(orderIndex).serviceCompletionTime = obj.clock;
            obj.orders(orderIndex).status = "staged";
            
            % Update service time statistic.
            if ~isnan(obj.orders(orderIndex).serviceStartTime)
                serviceDuration = obj.clock - obj.orders(orderIndex).serviceStartTime;
                serviceField = obj.statNames.warehouseServiceTime{warehouseId};
                obj.stats.update(serviceField, serviceDuration);
            end

            % Move order to the vehicle buffer.
            buffer = obj.state.vehicleBuffers{warehouseId};
            buffer(end+1) = orderIndex; 
            obj.state.vehicleBuffers{warehouseId} = buffer;

            % Start the next order if the queue is not empty.
            queue = obj.state.warehouseQueues{warehouseId};
            if ~isempty(queue)
                queueField = obj.statNames.warehouseQueueTimeAvg{warehouseId};
                oldQueueLength = obj.state.warehouseQueueLengths(warehouseId);
                obj.stats.update(queueField, oldQueueLength);
                nextOrder = queue(1);
                queue(1) = [];
                obj.state.warehouseQueues{warehouseId} = queue;
                obj.state.warehouseQueueLengths(warehouseId) = oldQueueLength - 1;
                obj.startWarehouseService(warehouseId, nextOrder);
            end

            % Check whether a capacity-triggered departure should be scheduled.
            obj.maybeScheduleVehicleDeparture(warehouseId);

        end %end handleWarehouseServiceCompletion

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function handleVehicleDeparture(obj, warehouseId, triggerType)
            %HANDLEVEHICLEDEPARTURE Dispatch a vehicle if possible.

            if nargin < 3 || isempty(triggerType)
                triggerType = 'schedule';
            end

            if strcmp(triggerType, 'capacity')
                obj.state.vehicleCapacityTriggerScheduled(warehouseId) = false;
            elseif strcmp(triggerType, 'schedule')
                obj.scheduleNextDepartureTimer(warehouseId, obj.clock);
            end

            if obj.isFleetFullyBusy(warehouseId)
                return; % All vehicles already on route.
            end

            buffer = obj.state.vehicleBuffers{warehouseId};
            if isempty(buffer)
                return;
            end

            capacity = obj.config.warehouses(warehouseId).vehicle.capacity;
            loadCount = min(capacity, numel(buffer));
            orderBatch = buffer(1:loadCount);
            obj.state.vehicleBuffers{warehouseId} = buffer(loadCount+1:end);

            % Update vehicle load factor statistic.
            loadField = obj.statNames.vehicleLoadFactor{warehouseId};
            obj.stats.update(loadField, loadCount);

            % Update vehicle utilisation (state change from idle to busy).
            utilField = obj.statNames.vehicleUtilization{warehouseId};
            obj.stats.update(utilField, obj.getVehicleBusyRatio(warehouseId));
            vehicleId = obj.acquireVehicleSlot(warehouseId);
            if vehicleId == 0
                % No idle vehicles available despite previous check; requeue the batch.
                obj.state.vehicleBuffers{warehouseId} = [orderBatch, obj.state.vehicleBuffers{warehouseId}];
                obj.maybeScheduleVehicleDeparture(warehouseId);
                return;
            end
            busyVector = obj.state.vehicleBusy{warehouseId};
            busyVector(vehicleId) = true;
            obj.state.vehicleBusy{warehouseId} = busyVector;

            % Update order status to reflect dispatch.
            for k = 1:numel(orderBatch)
                idx = orderBatch(k);
                obj.orders(idx).dispatchTime = obj.clock;
                obj.orders(idx).deliveryAttempts = obj.orders(idx).deliveryAttempts + 1;
                obj.orders(idx).status = "inTransit";
            end

            trip = obj.buildVehicleTrip(warehouseId, orderBatch);
            trip.vehicleId = vehicleId;
            fleetTrips = obj.state.activeTrips{warehouseId};
            fleetTrips{vehicleId} = trip;
            obj.state.activeTrips{warehouseId} = fleetTrips;

            % Schedule delivery events for each stop.
            for s = 1:numel(trip.stops)
                stop = trip.stops(s);
                obj.eventsList.Enqueue(LockerDeliveryEvent(warehouseId, vehicleId, s, stop.arrivalTime));
            end

            % Schedule the vehicle return event.
            obj.eventsList.Enqueue(VehicleReturnEvent(warehouseId, vehicleId, trip.returnTime));

            % After dispatching, check whether we can prepare another departure.
            obj.maybeScheduleVehicleDeparture(warehouseId);
        
        end %end handleVehicleDeparture
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function handleLockerDelivery(obj, warehouseId, vehicleId, stopIndex)
            %HANDLELOCKERDELIVERY Deliver a batch of parcels to a locker.

            fleetTrips = obj.state.activeTrips{warehouseId};
            if vehicleId < 1 || vehicleId > numel(fleetTrips)
                return;
            end

            trip = fleetTrips{vehicleId};

            if isempty(trip) || stopIndex > numel(trip.stops)
                return;
            end

            stop = trip.stops(stopIndex);
            lockerId = stop.lockerId;
            orderList = stop.orderIndices;
            if isempty(orderList)
                return;
            end

            occField = obj.statNames.lockerOccupancy{lockerId};
            oldOcc = obj.state.lockerOccupancy(lockerId);
            obj.stats.update(occField, oldOcc);

            capacity = obj.config.lockers(lockerId).capacity;
            available = capacity - oldOcc;
            toDeliver = numel(orderList);
            deliveredCount = min(toDeliver, max(0, available));
            blockedCount = toDeliver - deliveredCount;

            % Update blocking statistics per parcel attempt.
            blockField = obj.statNames.lockerBlocking{lockerId};
            for i = 1:deliveredCount
                obj.stats.update(blockField, 0);
            end
            for i = 1:blockedCount
                obj.stats.update(blockField, 1);
            end

            if deliveredCount > 0
                deliveredOrders = orderList(1:deliveredCount);
                newOcc = oldOcc + deliveredCount;
                obj.state.lockerOccupancy(lockerId) = newOcc;

                for k = 1:numel(deliveredOrders)
                    orderIdx = deliveredOrders(k);
                    obj.orders(orderIdx).deliveryTime = obj.clock;
                    obj.orders(orderIdx).status = "delivered";
                    leadTime = obj.clock - obj.orders(orderIdx).arrivalTime;
                    obj.stats.update(obj.statNames.leadTimeAverage, leadTime);
                    if obj.options.storeLeadTimes
                        if obj.leadTimeCount >= numel(obj.leadTimeSamples)
                            obj.leadTimeSamples(end+1, 1) = NaN;
                        end
                        obj.leadTimeCount = obj.leadTimeCount + 1;
                        obj.leadTimeSamples(obj.leadTimeCount) = leadTime;
                    end
                    obj.completedOrders = obj.completedOrders + 1;
                end

                % Track delivered orders inside the locker and schedule their
                % customer pickups using the order-specific delay when
                % available.
                lockerContents = obj.state.lockerContents{lockerId};
                lockerContents = [lockerContents, deliveredOrders(:).'];
                obj.state.lockerContents{lockerId} = lockerContents;

                for k = 1:numel(deliveredOrders)
                    orderIdx = deliveredOrders(k);
                    obj.scheduleLockerPickup(lockerId, orderIdx, obj.clock);
                end
            end

            if blockedCount > 0
                blockedOrders = orderList(deliveredCount+1:end);
                trip.undeliveredOrders = [trip.undeliveredOrders, blockedOrders(:).'];
                for k = 1:numel(blockedOrders)
                    orderIdx = blockedOrders(k);
                    obj.orders(orderIdx).status = "blocked";
                end
                obj.blockedParcels = obj.blockedParcels + blockedCount;
                obj.stats.update(obj.statNames.totalBlockedParcels, blockedCount);
            end

            trip.stops(stopIndex).delivered = deliveredCount;
            trip.stops(stopIndex).blocked = blockedCount;
            fleetTrips{vehicleId} = trip;
            obj.state.activeTrips{warehouseId} = fleetTrips;
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function handleVehicleReturn(obj, warehouseId, vehicleId)
            %HANDLEVEHICLERETURN Complete a vehicle trip and process undelivered
            %parcels (if any).
            
            busyVector = obj.state.vehicleBusy{warehouseId};
            if vehicleId < 1 || vehicleId > numel(busyVector)
                obj.maybeScheduleVehicleDeparture(warehouseId);
                return;
            end

            utilField = obj.statNames.vehicleUtilization{warehouseId};
            obj.stats.update(utilField, obj.getVehicleBusyRatio(warehouseId));
            busyVector(vehicleId) = false;
            obj.state.vehicleBusy{warehouseId} = busyVector;

            fleetTrips = obj.state.activeTrips{warehouseId};
            trip = fleetTrips{vehicleId};
            if isempty(trip)
                obj.maybeScheduleVehicleDeparture(warehouseId);
                return;
            end

            if ~isempty(trip.undeliveredOrders)
                buffer = obj.state.vehicleBuffers{warehouseId};
                % Reinsert blocked parcels at the front of the buffer to retry soon.
                obj.state.vehicleBuffers{warehouseId} = [trip.undeliveredOrders(:).', buffer];
                for k = 1:numel(trip.undeliveredOrders)
                    idx = trip.undeliveredOrders(k);
                    obj.orders(idx).status = "staged";
                end
            end

            fleetTrips{vehicleId} = [];
            obj.state.activeTrips{warehouseId} = fleetTrips;

            % After returning, check whether another departure is possible.
            obj.maybeScheduleVehicleDeparture(warehouseId);

        end  %end handleVehicleReturn

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function handleLockerPickup(obj, lockerId, orderIndex)
            %HANDLELOCKERPICKUP Manage the completion of a customer pickup.

            if lockerId < 1 || lockerId > obj.config.numLockers
                return;
            end

            occField = obj.statNames.lockerOccupancy{lockerId};
            oldOcc = obj.state.lockerOccupancy(lockerId);
            obj.stats.update(occField, oldOcc);

            if oldOcc <= 0
                obj.state.lockerOccupancy(lockerId) = 0;
                return;
            end

            contents = obj.state.lockerContents{lockerId};
            if isempty(contents)
                obj.state.lockerOccupancy(lockerId) = max(0, oldOcc - 1);
                return;
            end

            if nargin < 3
                orderIndex = 0;
            end

            if orderIndex > 0
                matchPos = find(contents == orderIndex, 1, 'first');
            else
                matchPos = find(contents <= 0, 1, 'first');
            end

            if isempty(matchPos)
                matchPos = 1;
            end

            contents(matchPos) = [];
            obj.state.lockerContents{lockerId} = contents;

            newOcc = max(0, oldOcc - 1);
            obj.state.lockerOccupancy(lockerId) = newOcc;

            if orderIndex > 0 && orderIndex <= numel(obj.orders)
                obj.orders(orderIndex).status = "pickedUp";
                obj.orders(orderIndex).pickupTime = obj.clock;
            end
        end

    end %end methods

    methods (Access = protected)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function startWarehouseService(obj, warehouseId, orderIndex)
            %STARTWAREHOUSESERVICE Begin service for the specified order.
            
            obj.state.warehouseInService(warehouseId) = 1;
            
            % Increase the number of orders currently in service at the
            % warehouse. The order identifier is already carried by the
            % scheduled completion event, therefore we just track the count.
            obj.state.warehouseInService(warehouseId) = obj.state.warehouseInService(warehouseId) + 1;

            obj.orders(orderIndex).serviceStartTime = obj.clock;
            obj.orders(orderIndex).status = "inService";

            waitField = obj.statNames.warehouseWaitTime{warehouseId};
            waitTime = obj.clock - obj.orders(orderIndex).arrivalTime;
            obj.stats.update(waitField, waitTime);

            serviceTime = obj.sampleDistribution(obj.config.warehouses(warehouseId).serviceTime);
            serviceTime = max(serviceTime, 0);
            completionTime = obj.clock + serviceTime;

            obj.eventsList.Enqueue(WarehouseServiceCompletionEvent(warehouseId, orderIndex, completionTime));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function maybeScheduleVehicleDeparture(obj, warehouseId)
            %MAYBESCHEDULEVEHICLEDEPARTURE Schedule a capacity-triggered departure
            %when the buffer is full and the vehicle is available.

            if obj.isFleetFullyBusy(warehouseId)
                return;
            end

            buffer = obj.state.vehicleBuffers{warehouseId};
            capacity = obj.config.warehouses(warehouseId).vehicle.capacity;
            if ~isfinite(capacity) || capacity <= 0
                capacityThreshold = 1;
            else
                capacityThreshold = capacity;
            end
            
            if numel(buffer) >= capacityThreshold && ~obj.state.vehicleCapacityTriggerScheduled(warehouseId)
                obj.state.vehicleCapacityTriggerScheduled(warehouseId) = true;
                obj.eventsList.Enqueue(VehicleDepartureEvent(warehouseId, 'capacity', obj.clock));
            end
        end %end maybeScheduleVehicleDeparture

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function vehicleId = acquireVehicleSlot(obj, warehouseId)
            %ACQUIREVEHICLESLOT Return the index of the first idle vehicle.
            busyVector = obj.state.vehicleBusy{warehouseId};
            vehicleId = find(~busyVector, 1, 'first');
            if isempty(vehicleId)
                vehicleId = 0;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function flag = isFleetFullyBusy(obj, warehouseId)
            %ISFLEETFULLYBUSY Check whether all vehicles from a warehouse are busy.
            busyVector = obj.state.vehicleBusy{warehouseId};
            if isempty(busyVector)
                flag = false;
            else
                flag = all(busyVector);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ratio = getVehicleBusyRatio(obj, warehouseId)
            %GETVEHICLEBUSYRATIO Fraction of vehicles currently busy.
            busyVector = obj.state.vehicleBusy{warehouseId};
            fleetSize = obj.config.warehouses(warehouseId).vehicle.fleetSize;
            if isempty(fleetSize) || fleetSize < 1
                fleetSize = 1;
            end
            denom = max(1, fleetSize);
            if isempty(busyVector)
                ratio = 0;
            else
                ratio = sum(busyVector) / denom;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function trip = buildVehicleTrip(obj, warehouseId, orderBatch)
            %BUILDVEHICLETRIP Create the route description for the dispatched batch.
            %   Vehicles follow a simple cyclic routing policy: each vehicle starts
            %   from its warehouse, visits lockers according to the predefined
            %   base route (skipping those without pending parcels) and finally
            %   returns to the warehouse. This lightweight policy enables several
            %   vehicles per depot to operate concurrently without solving a full
            %   routing optimisation problem.

            vehicleCfg = obj.config.warehouses(warehouseId).vehicle;
            assignedLockers = obj.config.warehouses(warehouseId).assignedLockers;
            if isempty(orderBatch)
                trip = struct('orders', orderBatch, 'stops', [], 'undeliveredOrders', [], 'returnTime', obj.clock);
                return;
            end
            if ~isfield(vehicleCfg,'originType') || isempty(vehicleCfg.originType)
                vehicleCfg.originType = 'warehouse';
            end
            if ~isfield(vehicleCfg,'returnToOrigin') || isempty(vehicleCfg.returnToOrigin)
                vehicleCfg.returnToOrigin = true;
            end
            if strcmpi(vehicleCfg.originType,'locker')
                if ~isfield(vehicleCfg,'originLockerId') || isempty(vehicleCfg.originLockerId)
                    error('originType is "locker" but originLockerId is not set.');
                end
                originLockerId = vehicleCfg.originLockerId;
            end
            batchSize = length(orderBatch);
            lockerIds = zeros(batchSize, 1);
            for i = 1:batchSize
                lockerIds(i) = obj.orders(orderBatch(i)).lockerId;
            end
            activeLockers = unique(lockerIds(:))';
            if isempty(activeLockers)
                trip = struct('orders', orderBatch, 'stops', [], 'undeliveredOrders', [], 'returnTime', obj.clock);
                return;
            end
            nStops = numel(activeLockers);
            N = nStops + 1;
            D = inf(N,N);
            for i = 1:N
                for j = 1:N
                    if i == j
                        continue
                    end
                    if i == 1 && j > 1
                        tgt = activeLockers(j-1);
                        if strcmpi(vehicleCfg.originType,'warehouse')
                            D(i,j) = max(obj.travelTimes.warehouseToLocker(warehouseId, tgt),0);
                        else
                            D(i,j) = max(obj.travelTimes.lockerToLocker(originLockerId, tgt),0);
                        end
                    elseif i > 1 && j == 1
                        src = activeLockers(i-1);
                        if strcmpi(vehicleCfg.originType,'warehouse')
                            if isfield(obj.travelTimes,'lockerToWarehouse')
                                D(i,j) = max(obj.travelTimes.lockerToWarehouse(src, warehouseId),0);
                            else
                                D(i,j) = max(obj.travelTimes.warehouseToLocker(warehouseId, src),0);
                            end
                        else
                            D(i,j) = max(obj.travelTimes.lockerToLocker(src, originLockerId),0);
                        end
                    else
                        src = activeLockers(i-1);
                        tgt = activeLockers(j-1);
                        D(i,j) = max(obj.travelTimes.lockerToLocker(src, tgt),0);
                    end
                end
            end
            idx = [];
            I = [];
            J = [];
            for i = 1:N
                for j = 1:N
                    if i ~= j
                        idx(end+1) = numel(idx)+1; 
                        I(end+1) = i;
                        J(end+1) = j;
                    end
                end
            end
            nArcs = numel(idx);
            f = zeros(nArcs,1);
            for k = 1:nArcs
                f(k) = D(I(k),J(k));
            end
            Aeq = zeros(2*N, nArcs);
            beq = ones(2*N,1);
            for i = 1:N
                row = i;
                cols = find(I == i);
                Aeq(row, cols) = 1;
            end
            for j = 1:N
                row = N + j;
                cols = find(J == j);
                Aeq(row, cols) = 1;
            end
            Aineq = [];
            bineq = [];
            lb = zeros(nArcs + (N-1),1);
            ub = ones(nArcs + (N-1),1);
            ub(nArcs+1:end) = N;
            lb(nArcs+1:end) = 2;
            intcon = 1:nArcs;
            for ii = 2:N
                for jj = 2:N
                    if ii == jj
                        continue
                    end
                    row = zeros(1, nArcs + (N-1));
                    arcCols = find(I == ii & J == jj);
                    if ~isempty(arcCols)
                        row(arcCols) = N;
                    end
                    row(nArcs + (ii-1)) = 1;
                    row(nArcs + (jj-1)) = -1;
                    Aineq = [Aineq; row];
                    bineq = [bineq; N-1]; 
                end
            end
            opt = optimoptions('intlinprog','Display','off');
            [sol,~,exitflag] = intlinprog(f,intcon, Aineq,bineq, Aeq,beq, lb,ub, opt);
            if exitflag <= 0
                baseRoute = assignedLockers(:)';
                orderPositions = arrayfun(@(x) find(baseRoute == x, 1, 'first'), activeLockers);
                for u = 1:numel(activeLockers)
                    if isempty(orderPositions(u))
                        baseRoute(end+1) = activeLockers(u);
                        orderPositions(u) = numel(baseRoute);
                    end
                end
                [~, sortIdx] = sort(orderPositions);
                stopOrder = activeLockers(sortIdx);
            else
                x = sol(1:nArcs) > 0.5;
                succ = zeros(1,N);
                for k = 1:nArcs
                    if x(k)
                        succ(I(k)) = J(k);
                    end
                end
                route = zeros(1,N+1);
                route(1) = 1;
                for t = 2:N+1
                    route(t) = succ(route(t-1));
                    if route(t) == 0
                        break
                    end
                    if route(t) == 1 && t < N+1
                        route = route(1:t);
                        break
                    end
                end
                tour = route(2:end-1);
                stopOrder = activeLockers(tour-1);
            end
            stops = repmat(struct('lockerId', 0, 'orderIndices', [], 'arrivalTime', 0, 'delivered', 0, 'blocked', 0), 1, numel(stopOrder));
            for s = 1:numel(stopOrder)
                locker = stopOrder(s);
                mask = (lockerIds == locker);
                stops(s).lockerId = locker;
                stops(s).orderIndices = orderBatch(mask);
            end
            currentTime = obj.clock;
            prevLocker = 0;
            for s = 1:numel(stops)
                locker = stops(s).lockerId;
                if prevLocker == 0
                    if strcmpi(vehicleCfg.originType,'warehouse')
                        travelTime = obj.travelTimes.warehouseToLocker(warehouseId, locker);
                    else
                        travelTime = obj.travelTimes.lockerToLocker(originLockerId, locker);
                    end
                else
                    travelTime = obj.travelTimes.lockerToLocker(prevLocker, locker);
                end
                currentTime = currentTime + max(travelTime, 0);
                stops(s).arrivalTime = currentTime;
                handlingBase = 0;
                if isfield(vehicleCfg, 'handlingTimePerStop') && ~isempty(vehicleCfg.handlingTimePerStop)
                    handlingBase = vehicleCfg.handlingTimePerStop;
                end
                handlingPerParcel = 0;
                if isfield(vehicleCfg, 'handlingTimePerParcel') && ~isempty(vehicleCfg.handlingTimePerParcel)
                    handlingPerParcel = vehicleCfg.handlingTimePerParcel;
                end
                lockerHandling = 0;
                if locker <= numel(obj.config.lockers) && isfield(obj.config.lockers(locker), 'handlingTime')
                    lockerHandling = obj.config.lockers(locker).handlingTime;
                end
                parcelsHere = numel(stops(s).orderIndices);
                serviceTime = handlingBase + handlingPerParcel * parcelsHere + lockerHandling;
                currentTime = currentTime + max(serviceTime, 0);
                prevLocker = locker;
            end
            if prevLocker == 0
                returnTime = currentTime;
            else
                if vehicleCfg.returnToOrigin
                    if strcmpi(vehicleCfg.originType,'warehouse')
                        if isfield(obj.travelTimes, 'lockerToWarehouse')
                            returnTravel = obj.travelTimes.lockerToWarehouse(prevLocker, warehouseId);
                        else
                            returnTravel = obj.travelTimes.warehouseToLocker(warehouseId, prevLocker);
                        end
                        returnTime = currentTime + max(returnTravel, 0);
                    else
                        returnTravel = obj.travelTimes.lockerToLocker(prevLocker, originLockerId);
                        returnTime = currentTime + max(returnTravel, 0);
                    end
                else
                    returnTime = currentTime;
                end
            end
            trip.orders = orderBatch;
            trip.stops = stops;
            trip.undeliveredOrders = [];
            trip.returnTime = returnTime;
        
        end %end buildVehicleTrip
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function scheduleNextDepartureTimer(obj, warehouseId, baseTime)
            %SCHEDULENEXTDEPARTURETIMER Schedule the next timetable departure.

            interval = obj.config.warehouses(warehouseId).vehicle.dispatchInterval;
            if isempty(interval) || ~isfinite(interval) || interval <= 0
                obj.state.vehicleNextScheduledDeparture(warehouseId) = Inf;
                return;
            end

            if nargin < 3
                baseTime = obj.clock;
            end

            nextTime = baseTime + interval;
            obj.state.vehicleNextScheduledDeparture(warehouseId) = nextTime;
            obj.eventsList.Enqueue(VehicleDepartureEvent(warehouseId, 'schedule', nextTime));
        
        end %end scheduleNextDepartureTimer

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function scheduleLockerPickup(obj, lockerId, orderIdx, currentTime)
            %SCHEDULELOCKERPICKUP Schedule the next customer pickup event.

            if nargin < 4
                currentTime = obj.clock;
            end
            
            pickupDelay = obj.orders(orderIdx).pickupDelay;
            nextTime = currentTime + pickupDelay * 60; 
            if ~isfinite(nextTime)
                return;
            end

            obj.eventsList.Enqueue(LockerPickupEvent(lockerId, orderIdx, nextTime));
        
        end %end scheduleLockerPickup
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = sampleDistribution(obj, distribution)
            %SAMPLEDISTRIBUTION Generic sampler for non-exponential distributions.

            if isempty(distribution)
                value = 0;
                return;
            end

            if isstruct(distribution)
                if isfield(distribution, 'sampleFcn') && ~isempty(distribution.sampleFcn)
                    if nargin(distribution.sampleFcn) == 1
                        value = distribution.sampleFcn(obj);
                    else
                        value = distribution.sampleFcn();
                    end
                    return;
                end
                if isfield(distribution, 'randomFcn') && ~isempty(distribution.randomFcn)
                    value = distribution.randomFcn();
                    return;
                end
                if isfield(distribution, 'value') && (~isfield(distribution, 'type') || any(strcmpi(distribution.type, {'deterministic', 'constant'})))
                    value = distribution.value;
                    return;
                end
                if ~isfield(distribution, 'type')
                    value = 0;
                    return;
                end
                type = lower(string(distribution.type));
            else
                value = distribution;
                return;
            end

            switch type
                case {'deterministic', 'constant'}
                    if isfield(distribution, 'value')
                        value = distribution.value;
                    elseif isfield(distribution, 'mean')
                        value = distribution.mean;
                    else
                        value = 0;
                    end
                case 'normal'
                    mu = 0;
                    sigma = 1;
                    if isfield(distribution, 'mu'); mu = distribution.mu; end
                    if isfield(distribution, 'mean'); mu = distribution.mean; end
                    if isfield(distribution, 'sigma'); sigma = distribution.sigma; end
                    if sigma < 0; sigma = abs(sigma); end
                    value = mu + sigma * obj.randnScalar();
                case 'lognormal'
                    mu = 0;
                    sigma = 1;
                    if isfield(distribution, 'mu'); mu = distribution.mu; end
                    if isfield(distribution, 'sigma'); sigma = distribution.sigma; end
                    if sigma < 0; sigma = abs(sigma); end
                    value = exp(mu + sigma * obj.randnScalar());
                case 'weibull'
                    shape = 1;
                    scale = 1;
                    if isfield(distribution, 'shape'); shape = distribution.shape; end
                    if isfield(distribution, 'scale'); scale = distribution.scale; end
                    u = max(min(obj.randUniform(), 1 - eps), eps);
                    value = scale * (-log(u))^(1/shape);
                case 'uniform'
                    a = 0; b = 1;
                    if isfield(distribution, 'a'); a = distribution.a; end
                    if isfield(distribution, 'b'); b = distribution.b; end
                    value = a + (b - a) * obj.randUniform();
                case 'exponential'
                    if isfield(distribution, 'rate') && distribution.rate > 0
                        value = obj.sampleExponential(distribution.rate);
                    elseif isfield(distribution, 'mean') && distribution.mean > 0
                        value = obj.sampleExponential(1/distribution.mean);
                    else
                        value = 0;
                    end
                case 'empirical'
                    if isfield(distribution, 'samples') && ~isempty(distribution.samples)
                        samples = distribution.samples(:);
                        idx = obj.randiScalar(numel(samples));
                        value = samples(idx);
                    else
                        value = 0;
                    end
                otherwise
                    error('Unsupported distribution type "%s".', distribution.type);
            end
        end %end sampleDistribution
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = sampleExponential(obj, rate)
            %SAMPLEEXPONENTIAL Sample exponential RV with the provided rate.
            u = max(min(obj.randUniform(), 1 - eps), eps);
            value = -log(u) / rate;
        end %end sampleExponential

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = randUniform(obj)
            %RANDUNIFORM Uniform random number in (0,1).
            if isempty(obj.randomStream)
                u = rand();
            else
                u = rand(obj.randomStream);
            end
        end %end randUniform

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z = randnScalar(obj)
            %RANDNSCALAR Standard normal random number.
            if isempty(obj.randomStream)
                z = randn();
            else
                z = randn(obj.randomStream);
            end
        end %end randnScalar

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function idx = randiScalar(obj, n)
            %RANDISCALAR Integer uniform random number between 1 and n.
            if isempty(obj.randomStream)
                idx = randi(n);
            else
                idx = randi(obj.randomStream, n);
            end
        end %end randiScalar

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function prepareStatisticNames(obj)
            %PREPARESTATISTICNAMES Pre-compute statistic identifiers.
            numW = obj.config.numWarehouses;
            numL = obj.config.numLockers;

            obj.statNames.warehouseQueueTimeAvg = cell(1, numW);
            obj.statNames.warehouseWaitTime = cell(1, numW);
            obj.statNames.warehouseServiceTime = cell(1, numW);
            obj.statNames.vehicleUtilization = cell(1, numW);
            obj.statNames.vehicleLoadFactor = cell(1, numW);
            for w = 1:numW
                obj.statNames.warehouseQueueTimeAvg{w} = sprintf('WarehouseQueueTimeAvg_%d', w);
                obj.statNames.warehouseWaitTime{w} = sprintf('WarehouseWaitTime_%d', w);
                obj.statNames.warehouseServiceTime{w} = sprintf('WarehouseServiceTime_%d', w);
                obj.statNames.vehicleUtilization{w} = sprintf('VehicleUtilization_%d', w);
                obj.statNames.vehicleLoadFactor{w} = sprintf('VehicleLoadFactor_%d', w);
            end

            obj.statNames.lockerOccupancy = cell(1, numL);
            obj.statNames.lockerBlocking = cell(1, numL);
            for l = 1:numL
                obj.statNames.lockerOccupancy{l} = sprintf('LockerOccupancy_%d', l);
                obj.statNames.lockerBlocking{l} = sprintf('LockerBlocking_%d', l);
            end

            obj.statNames.totalBlockedParcels = 'TotalBlockedParcels';
            obj.statNames.leadTimeAverage = 'LeadTimeAverage';
        
        end %end prepareStatisticNames
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function events = createInitialEvents(obj)
            %CREATEINITIALEVENTS Populate the initial future events list.
            events = event.empty;

            if obj.totalOrders >= 1
                events(end+1) = OrderArrivalEvent(1, obj.orders(1).arrivalTime);
            end

            for w = 1:obj.config.numWarehouses
                interval = obj.config.warehouses(w).vehicle.dispatchInterval;
                offset = [];
                if isfield(obj.config.warehouses(w).vehicle, 'initialOffset')
                    offset = obj.config.warehouses(w).vehicle.initialOffset;
                end
                if isempty(offset)
                    offset = interval;
                end
                if ~isempty(interval) && isfinite(interval) && interval > 0 && ~isempty(offset) && isfinite(offset) && offset > 0
                    obj.state.vehicleNextScheduledDeparture(w) = offset;
                    events(end+1) = VehicleDepartureEvent(w, 'schedule', offset);
                else
                    obj.state.vehicleNextScheduledDeparture(w) = Inf;
                end
            end

            for l = 1:obj.config.numLockers
                lockerContents = obj.state.lockerContents{l};
                if isempty(lockerContents)
                    continue;
                end
                baseTime = obj.clock;
                for idx = 1:numel(lockerContents)
                    orderIndex = lockerContents(idx);
                    nextTime = obj.computeNextPickupTime(l, orderIndex, baseTime);
                    if ~isfinite(nextTime)
                        continue;
                    end
                    events(end+1) = LockerPickupEvent(l, orderIndex, nextTime); %#ok<AGROW>
                    baseTime = max(baseTime, nextTime);
                end
            end
        end

    end %end methods (Access = protected)
        

    methods (Static, Access = private)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function parsed = parseModelConfig(modelConfig)
            if ~isstruct(modelConfig)
                error('modelConfig must be a struct.');
            end

            if ~isfield(modelConfig, 'orders')
                error('modelConfig must contain an ''orders'' field.');
            end
            orders = LastMileDelivery.prepareOrders(modelConfig.orders);

            if ~isfield(modelConfig, 'lockers')
                error('modelConfig must contain a ''lockers'' field.');
            end
            lockers = LastMileDelivery.prepareLockers(modelConfig.lockers);

            if ~isfield(modelConfig, 'warehouses')
                error('modelConfig must contain a ''warehouses'' field.');
            end
            warehouses = LastMileDelivery.prepareWarehouses(modelConfig.warehouses, ...
                orders, lockers, modelConfig.serviceTimesFiles, modelConfig.wCfg);

            numWarehouses = numel(warehouses);
            numLockers = numel(lockers);

            if ~isfield(modelConfig, 'travelTimes')
                error('modelConfig must contain a ''travelTimes'' field.');
            end
            travelTimes = LastMileDelivery.prepareTravelTimes(modelConfig.travelTimes, numWarehouses, numLockers);

            options = LastMileDelivery.defaultOptions();
            if isfield(modelConfig, 'options') && ~isempty(modelConfig.options)
                opts = modelConfig.options;
                fn = fieldnames(opts);
                for i = 1:numel(fn)
                    options.(fn{i}) = opts.(fn{i});
                end
            end

            randomStream = [];
            if isfield(modelConfig, 'randomStream') && ~isempty(modelConfig.randomStream)
                randomStream = modelConfig.randomStream;
            elseif isfield(modelConfig, 'rngStream') && ~isempty(modelConfig.rngStream)
                randomStream = modelConfig.rngStream;
            end
            if ~isempty(randomStream) && ~isa(randomStream, 'RandStream')
                error('randomStream must be a RandStream object.');
            end

            config.warehouses = warehouses;
            config.lockers = lockers;
            config.numWarehouses = numWarehouses;
            config.numLockers = numLockers;
            config.lockersByWarehouse = cell(1, numWarehouses);
            for w = 1:numWarehouses
                config.lockersByWarehouse{w} = warehouses(w).assignedLockers(:)';
            end

            parsed.orders = orders;
            parsed.config = config;
            parsed.travelTimes = travelTimes;
            parsed.options = options;
            parsed.randomStream = randomStream;
        end %end parseModelConfig
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function orders = prepareOrders(input)
            if istable(input)
                T = input;
            elseif isstruct(input)
                T = struct2table(input);
            else
                error('Unsupported type for orders. Use a table or struct array.');
            end

            varsLower = lower(T.Properties.VariableNames);

            arrivalIdx = find(strcmp(varsLower, 'timestampminutes'), 1);
            if isempty(arrivalIdx)
                arrivalIdx = find(strcmp(varsLower, 'arrivaltime'), 1);
            end
            if isempty(arrivalIdx)
                error('Orders data must contain an arrival time column (TimestampMinutes or ArrivalTime).');
            end
            arrivalTime = double(T{:, arrivalIdx});

            warehouseIdx = find(strcmp(varsLower, 'warehouse'), 1);
            if isempty(warehouseIdx)
                error('Orders data must contain a Warehouse column.');
            end
            warehouseId = double(T{:, warehouseIdx});

            lockerIdx = find(strcmp(varsLower, 'locker_id'), 1);
            if isempty(lockerIdx)
                lockerIdx = find(strcmp(varsLower, 'locker'), 1);
            end
            if isempty(lockerIdx)
                error('Orders data must contain a Locker_ID column.');
            end
            lockerRaw = T{:, lockerIdx};
            if iscell(lockerRaw) || isstring(lockerRaw)
                lockerNum = zeros(size(lockerRaw));
                for i = 1:numel(lockerRaw)
                    token = regexprep(char(lockerRaw(i)), '\D', '');
                    lockerNum(i) = str2double(token);
                end
            else
                lockerNum = double(lockerRaw);
            end

            idIdx = find(strcmp(varsLower, 'order_id'), 1);
            if isempty(idIdx)
                idData = arrayfun(@(k) sprintf('ORD_%05d', k), (1:height(T))', 'UniformOutput', false);
            else
                raw = T{:, idIdx};
                idData = cellstr(string(raw));
            end

            pickupIdx = find(contains(varsLower, 'pickup_delay'), 1);
            if ~isempty(pickupIdx)
                pickupDelay = double(T{:, pickupIdx});
            else
                pickupDelay = nan(height(T), 1);
            end

            n = numel(arrivalTime);
            [arrivalTime, sortIdx] = sort(arrivalTime(:));
            warehouseId = warehouseId(sortIdx);
            lockerNum = lockerNum(sortIdx);
            idData = idData(sortIdx);
            pickupDelay = pickupDelay(sortIdx);

            orders = repmat(struct('id', '', 'arrivalTime', 0, 'warehouseId', 0, 'lockerId', 0, ...
                'pickupDelay', NaN, 'pickupTime', NaN, 'serviceStartTime', NaN, 'serviceCompletionTime', NaN, ...
                'dispatchTime', NaN, 'deliveryTime', NaN, 'deliveryAttempts', 0, 'status', "created"), n, 1);
            
            for i = 1:n
                orders(i).id = idData{i};
                orders(i).arrivalTime = arrivalTime(i);
                orders(i).warehouseId = warehouseId(i);
                orders(i).lockerId = lockerNum(i);
                orders(i).pickupDelay = pickupDelay(i);
            end
        end %end prepareOrders
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function warehouses = prepareWarehouses(input, orders, lockers, files, cfg)
            if istable(input)
                raw = table2struct(input);
            elseif isstruct(input)
                raw = input;
            else
                error('warehouses must be provided as a table or struct array.');
            end

            numW = numel(raw);
            warehouses = repmat(struct('id', '', 'serviceTime', [], 'vehicle', struct(), 'assignedLockers', [], 'name', ''), numW, 1);

            lockerIds = arrayfun(@(l) l.id, lockers);
            
            ex = readtable(files{1});
            new = readtable(files{2});

            for w = 1:numW
                entry = raw(w);
                if isfield(entry, 'id') && ~isempty(entry.id)
                    warehouses(w).id = char(string(entry.id));
                else
                    warehouses(w).id = sprintf('W%d', w);
                end
                if isfield(entry, 'name') && ~isempty(entry.name)
                    warehouses(w).name = char(string(entry.name));
                else
                    warehouses(w).name = warehouses(w).id;
                end

                if isfield(entry, 'serviceTime') && ~isempty(entry.serviceTime)
                    warehouses(w).serviceTime = entry.serviceTime;
                elseif isfield(entry, 'service_time') && ~isempty(entry.service_time)
                    warehouses(w).serviceTime = entry.service_time;
                elseif isfield(entry, 'serviceMean') && ~isempty(entry.serviceMean)
                    warehouses(w).serviceTime = struct('type', 'deterministic', 'value', entry.serviceMean);
                else
                    id = cfg(w);
                    if id < 3
                        warehouses(w).serviceTime = struct(...
                            'type', 'lognormal', ...
                            'mu', log(ex{id, 3} * 60), ...
                            'sigma', 0.5);
                    else
                        warehouses(w).serviceTime = struct(...
                            'type', 'lognormal', ...
                            'mu', log(new{id - 2, 5} * 60), ...
                            'sigma', 0.5);
                    end
                end

                if isfield(entry, 'vehicle') && ~isempty(entry.vehicle)
                    vehicle = entry.vehicle;
                else
                    vehicle = struct();
                end
                if ~isfield(vehicle, 'capacity') || isempty(vehicle.capacity)
                    error('Each warehouse must define vehicle.capacity.');
                end
                if ~isfield(vehicle, 'dispatchInterval')
                    vehicle.dispatchInterval = Inf;
                end
                if ~isfield(vehicle, 'handlingTimePerStop')
                    vehicle.handlingTimePerStop = 0;
                end
                if ~isfield(vehicle, 'handlingTimePerParcel')
                    vehicle.handlingTimePerParcel = 1/7;
                end
                if ~isfield(vehicle, 'initialOffset')
                    vehicle.initialOffset = vehicle.dispatchInterval;
                end
                if isfield(vehicle, 'fleetSize') && ~isempty(vehicle.fleetSize)
                    fleetSize = vehicle.fleetSize;
                elseif isfield(vehicle, 'numVehicles') && ~isempty(vehicle.numVehicles)
                    fleetSize = vehicle.numVehicles;
                elseif isfield(vehicle, 'vehicles') && ~isempty(vehicle.vehicles)
                    fleetSize = vehicle.vehicles;
                else
                    fleetSize = 1;
                end
                if ~isscalar(fleetSize) || isnan(fleetSize)
                    fleetSize = 1;
                end
                fleetSize = max(1, round(double(fleetSize)));
                vehicle.fleetSize = fleetSize;

                warehouses(w).vehicle = vehicle;

                if isfield(entry, 'assignedLockers') && ~isempty(entry.assignedLockers)
                    assigned = double(entry.assignedLockers(:)');
                else
                    mask = arrayfun(@(o) o.warehouseId == w, orders);
                    assigned = unique(arrayfun(@(o) o.lockerId, orders(mask)));
                end
                assigned = assigned(:)';
                assigned = assigned(ismember(assigned, lockerIds));
                warehouses(w).assignedLockers = assigned;

                if ~isfield(warehouses(w).vehicle, 'route') || isempty(warehouses(w).vehicle.route)
                    warehouses(w).vehicle.route = assigned;
                else
                    route = double(warehouses(w).vehicle.route(:)');
                    missing = setdiff(assigned, route, 'stable');
                    warehouses(w).vehicle.route = [route, missing];
                end
            end
        end %end prepareWarehouses

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lockers = prepareLockers(input)
            if istable(input)
                raw = table2struct(input);
            elseif isstruct(input)
                raw = input;
            else
                error('lockers must be provided as a table or struct array.');
            end

            numL = numel(raw);
            lockers = repmat(struct('id', 0, 'capacity', 0, 'initialOccupancy', 0, 'pickupRate', 1/60, ...
                'pickupServiceTime', struct('type', 'deterministic', 'value', 0), 'handlingTime', 0, 'name', ''), numL, 1);

            for l = 1:numL
                entry = raw(l);
                if isfield(entry, 'id') && ~isempty(entry.id)
                    lockers(l).id = double(entry.id);
                elseif isfield(entry, 'Locker_ID') && ~isempty(entry.Locker_ID)
                    s = string(entry.Locker_ID);
                    d = regexp(s, '\d+', 'match', 'once');
                    if ~isempty(d)
                        lockers(l).id = str2double(d);
                    else
                        lockers(l).id = l;             
                    end
                else
                    lockers(l).id = l;
                end

                if isfield(entry, 'name') && ~isempty(entry.name)
                    lockers(l).name = char(string(entry.name));
                else
                    lockers(l).name = sprintf('Locker_%d', lockers(l).id);
                end

                if isfield(entry, 'Capacity') && ~isempty(entry.Capacity)
                    lockers(l).capacity = double(entry.Capacity);
                elseif isfield(entry, 'capacity') && ~isempty(entry.capacity)
                    lockers(l).capacity = double(entry.capacity);
                else
                    error('Each locker must define a capacity.');
                end

                if isfield(entry, 'initialOccupancy') && ~isempty(entry.initialOccupancy)
                    lockers(l).initialOccupancy = min(double(entry.initialOccupancy), lockers(l).capacity);
                else
                    lockers(l).initialOccupancy = 0;
                end

                if isfield(entry, 'pickupRate') && ~isempty(entry.pickupRate)
                    lockers(l).pickupRate = double(entry.pickupRate);
                elseif isfield(entry, 'meanPickupInterarrival') && ~isempty(entry.meanPickupInterarrival)
                    lockers(l).pickupRate = 1 / double(entry.meanPickupInterarrival);
                else
                    lockers(l).pickupRate = 1 / (5 * 60); % default: one pickup every 12 hours
                end

                if isfield(entry, 'pickupServiceTime') && ~isempty(entry.pickupServiceTime)
                    lockers(l).pickupServiceTime = entry.pickupServiceTime;
                elseif isfield(entry, 'pickupService') && ~isempty(entry.pickupService)
                    lockers(l).pickupServiceTime = entry.pickupService;
                else
                    lockers(l).pickupServiceTime = struct('type', 'deterministic', 'value', 0);
                end

                if isfield(entry, 'handlingTime') && ~isempty(entry.handlingTime)
                    lockers(l).handlingTime = double(entry.handlingTime);
                else
                    lockers(l).handlingTime = 0;
                end
            end
        end %end prepareLockers

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function travel = prepareTravelTimes(input, numWarehouses, numLockers)
            if ~isstruct(input)
                error('travelTimes must be a struct with matrix fields.');
            end

            if ~isfield(input, 'warehouseToLocker')
                error('travelTimes must contain a warehouseToLocker field.');
            end
            w2l = LastMileDelivery.toNumericMatrix(input.warehouseToLocker);
            if size(w2l, 1) == numLockers && size(w2l, 2) == numWarehouses
                w2l = w2l';
            end
            if size(w2l, 1) ~= numWarehouses || size(w2l, 2) ~= numLockers
                error('warehouseToLocker matrix has inconsistent dimensions.');
            end

            if isfield(input, 'lockerToLocker')
                l2l = LastMileDelivery.toNumericMatrix(input.lockerToLocker);
                if size(l2l, 1) ~= numLockers || size(l2l, 2) ~= numLockers
                    error('lockerToLocker matrix has inconsistent dimensions.');
                end
            else
                l2l = zeros(numLockers);
            end

            if isfield(input, 'lockerToWarehouse')
                l2w = LastMileDelivery.toNumericMatrix(input.lockerToWarehouse);
                if size(l2w, 1) == numWarehouses && size(l2w, 2) == numLockers
                    l2w = l2w';
                end
                if size(l2w, 1) ~= numLockers || size(l2w, 2) ~= numWarehouses
                    error('lockerToWarehouse matrix has inconsistent dimensions.');
                end
            else
                l2w = w2l';
            end

            travel.warehouseToLocker = w2l;
            travel.lockerToLocker = l2l;
            travel.lockerToWarehouse = l2w;
        
        end %end prepareTravelTimes
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function stateStruct = createInitialState(parsed)
            numWarehouses = parsed.config.numWarehouses;
            numLockers = parsed.config.numLockers;

            warehouseQueues = cell(1, numWarehouses);
            warehouseQueueLengths = zeros(1, numWarehouses);
            warehouseInService = zeros(1, numWarehouses);
            vehicleBuffers = cell(1, numWarehouses);
            vehicleBusy = cell(1, numWarehouses);
            vehicleCapacityTriggerScheduled = false(1, numWarehouses);
            vehicleNextScheduledDeparture = inf(1, numWarehouses);
            activeTrips = cell(1, numWarehouses);
            lockerOccupancy = zeros(1, numLockers);
            lockerContents = cell(1, numLockers);

            for w = 1:numWarehouses
                warehouseQueues{w} = [];
                vehicleBuffers{w} = [];
                fleetSize = parsed.config.warehouses(w).vehicle.fleetSize;
                if isempty(fleetSize) || fleetSize < 1
                    fleetSize = 1;
                end
                vehicleBusy{w} = false(1, fleetSize);
                activeTrips{w} = cell(1, fleetSize);
            end
            for l = 1:numLockers
                initialOcc = parsed.config.lockers(l).initialOccupancy;
                lockerOccupancy(l) = initialOcc;
                if initialOcc > 0
                    lockerContents{l} = zeros(1, initialOcc);
                else
                    lockerContents{l} = [];
                end
            end

            stateStruct.warehouseQueues = warehouseQueues;
            stateStruct.warehouseQueueLengths = warehouseQueueLengths;
            stateStruct.warehouseInService = warehouseInService;
            stateStruct.vehicleBuffers = vehicleBuffers;
            stateStruct.vehicleBusy = vehicleBusy;
            stateStruct.vehicleCapacityTriggerScheduled = vehicleCapacityTriggerScheduled;
            stateStruct.vehicleNextScheduledDeparture = vehicleNextScheduledDeparture;
            stateStruct.activeTrips = activeTrips;
            stateStruct.lockerOccupancy = lockerOccupancy;
            stateStruct.lockerContents = lockerContents;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function statStruct = createStatistics(parsed)
            numWarehouses = parsed.config.numWarehouses;
            numLockers = parsed.config.numLockers;

            statStruct = struct();
            for w = 1:numWarehouses
                statStruct.(sprintf('WarehouseQueueTimeAvg_%d', w)) = TimeAverageStatistic();
                statStruct.(sprintf('WarehouseWaitTime_%d', w)) = AverageStatistic();
                statStruct.(sprintf('WarehouseServiceTime_%d', w)) = AverageStatistic();
                statStruct.(sprintf('VehicleUtilization_%d', w)) = TimeAverageStatistic();
                statStruct.(sprintf('VehicleLoadFactor_%d', w)) = AverageStatistic();
            end
            for l = 1:numLockers
                statStruct.(sprintf('LockerOccupancy_%d', l)) = TimeAverageStatistic();
                statStruct.(sprintf('LockerBlocking_%d', l)) = AverageStatistic();
            end
            statStruct.TotalBlockedParcels = SumStatistic();
            statStruct.LeadTimeAverage = AverageStatistic();
        end %end createStatistics

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function options = defaultOptions()
            options.storeLeadTimes = true;
            options.stopWhenAllDelivered = true;
        end %end defaultOptions
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function matrix = toNumericMatrix(data)
            if istable(data)
                matrix = table2array(data);
            elseif isnumeric(data)
                matrix = double(data);
            elseif iscell(data)
                matrix = cellfun(@double, data);
            else
                error('Unsupported matrix representation.');
            end
        end %end toNumericMatrix

    end % end methods (Static, Access = private)

end %end classdef