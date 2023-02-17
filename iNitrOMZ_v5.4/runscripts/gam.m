% HERE, specify iNitrOMZ root path ($PATHTOINSTALL/iNitrOMZ/)
bgc1d_root='/Users/yangsi/Box Sync/UCLA/MATLAB/BGC_Bianchi/iNitrOMZ/';
addpath(genpath(bgc1d_root)); % adds root to MATLAB's search paths

% Parameters to tune.
remin = 0.08/86400;
param_name1 = {'wup_param'   'poc_flux_top'    'KO2Rem'   'KNO2No' 'KO2Den1'  'KO2Den2' };
param_min1  = [  0.5e-7         -30/86400       0.01       0.01       0.01      0.01  ];
param_max1  = [ 15.0e-7          -3/86400        3.0      0.5        8.0        8.0  ];

param_name2 = {'KO2Den3'   'KNO3Den1'   'KDen1'  'KDen2'   'KDen3'};
param_min2  = [   0.01       0.01      remin/10  remin/10  remin/10];
param_max2  = [    8          8         remin      remin    remin ];

param_name3 = {'KNO2Den2' 'KN2ODen3'  'KAx'     'KNH4Ax'  'KNO2Ax'  'KO2Ax'};
param_min3  = [  0.001      0.001    remin/10     0.01      0.01      0.01 ];
param_max3  = [   0.5        0.5      remin*2     0.5      0.5      8  ];

% Constraints: KDen1 + KDen2 + KDen3 = remin
Aeq = [0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0];
beq = remin;

% concatenate names, mins, maxs
param = [param_name1 param_name2 param_name3];
param_min = [param_min1 param_min2 param_min3];
param_max = [param_max1 param_max2 param_max3];
nparam = size(param,2); % number of parameters

% make handle for cost function. x is an array of parameter value, param are the names which need to be passed 
costfunc = @(x)bgc1d_fc2minimize(x,param);

% Options
options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf,'UseParallel', true, 'UseVectorized', false,'Display','iter');

% Optimization
x = ga(costfunc,nparam,[],[],Aeq,beq,param_min,param_max,[],options);

% Save ga output using today's date
DateNow = bgc1d_getDate();
save([bgc1d_root,'saveOut/ga_out_',DateNow,'.mat']);

