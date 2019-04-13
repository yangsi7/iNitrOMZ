% ========================================================================
%
% File     : Cost.m
% Date     : September 2017
% Author   : Simon Yang
% Function : compute individual cost function 
% ========================================================================
% Quadratic cost, normalized at each depth by the standard deviation of the data, overall cost normalized by number of ponts
% sum((z-z_data)^2/std^2)/n_data


function Cost = Cost_quad_std_weighted(scenario_in, Data)
    % If no standard deviation, use 20% of value at that point
    [pred, dat] = minmax(scenario_in, Data);
    Cost = nansum((pred - dat.val).^2,2)./nansum(~isnan(dat.val),2);

