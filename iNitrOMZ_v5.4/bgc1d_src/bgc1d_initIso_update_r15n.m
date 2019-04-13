function bgc = bgc1d_initIso_update_r15n(bgc,t)

 % % % % % % % % % % % % % % % % % % % % % % % % % %
 % This function updates isotope ratios 
 % % % % % % % % % % % % % % % % % % % % % % % % % %

 %%%%%  ratios  %%%%%
 %%%%%%%%%%%%%%%%%%%%

 %% Initialize isotopic ratios %%
 bgc.r15norg = ones(size(bgc.zgrid));
 bgc.r15no3 = ones(size(bgc.zgrid));
 bgc.r15no2 = ones(size(bgc.zgrid));
 bgc.r15nh4 = ones(size(bgc.zgrid));
 bgc.r15n2o = ones(size(bgc.zgrid));

% %%%%%%%% Corg %%%%%%%
% idx = t.poc>0;
% ii = dNiso('d15N', bgc.d15norg, 'N', t.poc.*bgc.NCrem);
% bgc.r15norg(idx) = ii.i15N(idx)./ii.N(idx);     % 15N/(14N+15N) of organic nitrogen
% bgc.r15norg(~idx) = 1.0;
%
% %%%%%%%% NO3 %%%%%%%%
% idx = t.no3>0;
% bgc.r15no3(idx) = t.i15no3(idx)./t.no3(idx);     % 15N/(14N+15N) of nitrate
% bgc.r15no3(~idx) = 1.0;
%
% %%%%%%%% NO2 %%%%%%%%
% idx = t.no2>0;
% bgc.r15no2(idx) = t.i15no2(idx)./t.no2(idx);     % 15N/(14N+15N) of nitrite
% bgc.r15no2(~idx) = 1.0;
%
% %%%%%%%% NH4 %%%%%%%%
% idx = t.nh4>0;
% bgc.r15nh4(idx) = t.i15nh4(idx)./t.nh4(idx);     % 15N/(14N+15N) of ammonia
% bgc.r15nh4(~idx) = 1.0;
%
% %%%%%%%% N2O %%%%%%%%
% idx = t.n2o>0;
% i15n2o = t.i15n2oA + t.i15n2oB;
% bgc.r15n2o(idx) = i15n2o(idx)./t.n2o(idx);     % 15N/(14N+15N) of nitrous oxide
% bgc.r15n2o(~idx) = 1.0;
eps=10^-23;

 %%%%%%%% Corg %%%%%%%
 ii = dNiso('d15N', bgc.d15norg, 'N', t.poc.*bgc.NCrem);
 bgc.r15norg = (ii.i15N+eps)./(ii.N+eps);     % 15N/(14N+15N) of organic nitrogen

 %%%%%%%% NO3 %%%%%%%%
 bgc.r15no3 = (t.i15no3+eps)./(t.no3+eps);     % 15N/(14N+15N) of nitrate

 %%%%%%%% NO2 %%%%%%%%
 bgc.r15no2 = (t.i15no2+eps)./(t.no2+eps);     % 15N/(14N+15N) of nitrite

 %%%%%%%% NH4 %%%%%%%%
 bgc.r15nh4 = (t.i15nh4+eps)./(t.nh4+eps);     % 15N/(14N+15N) of ammonia

 %%%%%%%% N2O %%%%%%%%
 i15n2o = t.i15n2oA + t.i15n2oB;
 bgc.r15n2o = (i15n2o+eps)./(t.n2o+eps);     % 15N/(14N+15N) of nitrous oxide
