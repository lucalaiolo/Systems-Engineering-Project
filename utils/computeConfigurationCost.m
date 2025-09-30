function totalCost = computeConfigurationCost(cfg, existingData, newData)
%COMPUTECONFIGURATIONCOST Compute total cost (capex + opex) for a configuration.
%   cfg is a column vector containing the warehouse indices that form the
%   configuration. existingData and newData are numeric matrices whose rows
%   correspond to existing and new warehouses respectively. The function
%   assumes that existing warehouses have no additional capex requirements
%   and that their opex is stored in column 12 of existingData. For new
%   warehouses, the opex is expected in column 15 and the capex in column 16
%   of newData.

numExisting = size(existingData, 1);
totalCost = 0;

for wIdx = reshape(cfg, 1, [])
    if wIdx <= numExisting
        operatingExpenditure = existingData(wIdx, 12);
        capitalExpenditure = 0;
    else
        newIdx = wIdx - numExisting;
        operatingExpenditure = newData(newIdx, 15);
        capitalExpenditure = newData(newIdx, 16);
    end

    totalCost = totalCost + operatingExpenditure + capitalExpenditure;
end

end