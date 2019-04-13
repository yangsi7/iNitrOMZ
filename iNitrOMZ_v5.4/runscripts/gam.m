cd('/oasis/scratch/comet/yangsi/temp_project/MATLAB/NITRO_MGA_v5.2/src');
remin = 0.08/86400;
param_name1 = {'KO2Rem'   'KNO2No' 'KO2Den1'  'KO2Den2' };
param_min1  = [ 0.01       0.01       0.01      0.01  ];
param_max1  = [   3.0      0.5        8.0        8.0  ];

param_name2 = {'KO2Den3'   'KNO3Den1'   'KDen1'  'KDen2'   'KDen3'};
param_min2  = [   0.01       0.01      remin/10  remin/10  remin/10];
param_max2  = [    8          8         remin      remin    remin ];

param_name3 = {'KNO2Den2' 'KN2ODen3'  'KAx'     'KNH4Ax'  'KNO2Ax'  'KO2Ax'};
param_min3  = [  0.001      0.001    remin/10     0.01      0.01      0.01   ];
param_max3  = [   0.5        0.5      remin*2     0.5      0.5      8  ];

Aeq = [0 0 0 0 0 0 1 1 1 0 0 0 0 0 0];
beq = remin;

param = [param_name1 param_name2 param_name3];
param_min = [param_min1 param_min2 param_min3];
param_max = [param_max1 param_max2 param_max3];

nparam = size(param,2);

options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf,'UseParallel', true, 'UseVectorized', false,'Display','iter');

x = ga(@bgc1d_fc2minimize,nparam,[],[],Aeq,beq,param_min,param_max,[],options);

