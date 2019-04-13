bgc.root='/Users/yangsi/Box Sync/UCLA/MATLAB/BGC_Bianchi/iNitrOMZ/iNitrOMZ_v5.4/';
cd(bgc.root);
bgc1d_paths_init;

tic;
bgc = bgc1d_initialize();
Tracer.name = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2'};
if strcmp(bgc.region,'ETNP')
	load([bgc.root,'/Data/compilation_ETNP_gridded.mat']);
	Data.names = Tracer.name;
	Data = GA_data_init(bgc,compilation_ETNP_gridded,Tracer.name);
elseif strcmp(bgc.region,'ETSP')
	load([bgc.root,'/Data/compilation_ETSP_gridded_Feb232018.mat']);
	Data.names = Tracer.name;
	Data = GA_data_init_opt(bgc,compilation_ETSP_gridded,Tracer.name);
end

if bgc.flux_diag == 1 
	[bgc.sol_time bgc.adv_time bgc.diff_time bgc.sms_time bgc.rest_time] = bgc1d_advection_diff(bgc);
	bgc.sol = squeeze(bgc.sol_time(end,:,:));
	bgc.sadv = squeeze(bgc.adv_time(end,:,:));
	bgc.sdiff = squeeze(bgc.diff_time(end,:,:));
	bgc.ssms = squeeze(bgc.sms_time(end,:,:));
	bgc.srest = squeeze(bgc.rest_time(end,:,:));
else
	[bgc.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc);
end

toc;
bgc = bgc1d_postprocess(bgc, Data);
bgc1d_plot(bgc);
