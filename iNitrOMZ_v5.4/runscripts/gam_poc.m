cd('/oasis/scratch/comet/yangsi/temp_project/MATLAB/NITRO_MGA_v5.2/src');
param_min = [-25.0/86400 0.5*10^-7 -0.9];
param_max = [-3.0/86400 3.0 *10^-7 -0.3];

nparam = 3;

options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf,'UseParallel', true, 'UseVectorized', false,'Display','iter');

[x,Fval,exitFlag,Output] = ga(@bgc1d_fc2minimize_poc,nparam,[],[],[],[],param_min,param_max,[],options);

save('/oasis/scratch/comet/yangsi/temp_project/MATLAB/NITRO_MGA_v5.2/GA_output/gamarch4_poc.mat','x');
