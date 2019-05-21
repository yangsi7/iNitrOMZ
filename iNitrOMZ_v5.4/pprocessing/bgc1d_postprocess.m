 function bgc = bgc1d_postprocess(bgc, Data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d Ncycle v1.0 - Simon Yang  - April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the ini solution from the bgc.sol structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     zstep = bgc.npt;
     bgc.z = linspace(bgc.zbottom,bgc.ztop,zstep+1);
     bgc.dz = (bgc.ztop - bgc.zbottom) / zstep;
     
     bgc.SolNames = {'o2','no3','poc','po4','n2o', 'nh4', 'no2', 'n2'};
     if bgc.RunIsotopes
	     bgc.SolNames = {'o2','no3','poc','po4','n2o', 'nh4', ...
	     'no2', 'n2', 'i15no3', 'i15no2', 'i15nh4', 'i15n2oA', 'i15n2oB' };
     end
     bgc.sol = squeeze(bgc.sol_time(end,:,:));
     ntrSol = length(bgc.SolNames);
     for indt=1:ntrSol
     	bgc.(bgc.SolNames{indt}) = bgc.sol(indt,:);
     	if bgc.flux_diag == 1
     	        bgc.(['adv' bgc.SolNames{indt}]) = bgc.sadv(indt,:);
     	        bgc.(['diff' bgc.SolNames{indt}]) = bgc.sdiff(indt,:);
     	        bgc.(['sms' bgc.SolNames{indt}]) = bgc.ssms(indt,:);
     	        bgc.(['rest' bgc.SolNames{indt}]) = bgc.srest(indt,:);
     	end
     	bgc.(['d' bgc.SolNames{indt}]) = nan(size(bgc.(bgc.SolNames{indt})));
     	bgc.(['d' bgc.SolNames{indt}])(2:end-1) = ...
	(bgc.(bgc.SolNames{indt})(3:end) - bgc.(bgc.SolNames{indt})(1:end-2))/(-2*bgc.dz);
     	bgc.(['d' bgc.SolNames{indt}])(1) = 0; 
     	bgc.(['d' bgc.SolNames{indt}])(end) = 0;
     	bgc.(['d2' bgc.SolNames{indt}]) = nan(size(bgc.(['d' bgc.SolNames{indt}])));
     	bgc.(['d2' bgc.SolNames{indt}])(2:end-1) = ...
	(bgc.(bgc.SolNames{indt})(3:end) - 2 * bgc.(bgc.SolNames{indt})(2:end-1) ...
	+ bgc.(bgc.SolNames{indt})(1:end-2))/(bgc.dz^2);
     	bgc.(['d2' bgc.SolNames{indt}])(1) = 0;
     	bgc.(['d2' bgc.SolNames{indt}])(end) = 0;
    end
    if bgc.RunIsotopes
	    eps = 10^-3;
	    ii = dNiso('i15N', bgc.i15no3, 'N', bgc.no3);
	    idx = find(bgc.no3<eps | bgc.i15no3<0);
	    bgc.d15no3 = ii.d15N;
	    bgc.d15no3(idx)=nan;
	    ii = dNiso('i15N', bgc.i15no2, 'N', bgc.no2);
            idx = find(bgc.no2<eps | bgc.i15no2<0);
	    bgc.d15no2 = ii.d15N;
	    bgc.d15no2(idx)=nan;
	    ii = dNiso('i15N', bgc.i15nh4, 'N', bgc.nh4);
	    idx = find(bgc.nh4<eps | bgc.i15no2<0);
	    bgc.d15nh4 = ii.d15N;
	    bgc.d15nh4(idx)=nan;
	    ii = dNiso('i15N', bgc.i15n2oA, 'N', bgc.n2o);
	    idx = find(bgc.n2o<eps/1000 | bgc.i15n2oA<0);
	    bgc.d15n2oA = ii.d15N;
	    bgc.d15n2oA(idx)=nan;
	    ii = dNiso('i15N', bgc.i15n2oB, 'N', bgc.n2o);
	    idx = find(bgc.n2o<eps/1000 | bgc.i15n2oB<0);
	    bgc.d15n2oB = ii.d15N;
	    bgc.d15n2oB(idx)=nan;
    end

    ntrData = length(bgc.SolNames);
    tmp = strcat('Data_', bgc.SolNames);
    for indt=1:ntrData
	try 
    	    bgc.(tmp{indt}) = Data.val(indt,:);
	catch 
            display(['Warning: did not find', tmp{indt}]);
	end
    end
    
 % Additional tracers:
 % NSTAR - NO3 deficit versus PO4
 bgc.nstar = bgc.no3-bgc.NCrem/bgc.PCrem*bgc.po4;

 bgc.dwup = nan(size(bgc.wup)); bgc.dwup(1)=0; bgc.dwup(end)=0;
 bgc.dwup(2:end-1) = (bgc.wup(3:end) - bgc.wup(1:end-2))/(-2*bgc.dz);
 bgc.ssN2OAdv =  - (bgc.wup .* bgc.dn2o + bgc.dwup .* bgc.n2o);
 bgc.ssN2ODiff =bgc.Kv .* bgc.d2n2o; %bgc.vdKv .* bgc.dn2o + bgc.Kv .* bgc.d2n2o;
 bgc.sms_n2o(1)=0;bgc.sms_n2o(end)=0;bgc.n2o_rest(1) =0; bgc.n2o_rest(end)=0;
 % BGC sources and sink
 % Biological sources and sinks terms
 sms =  bgc1d_sms_diag(bgc); 

 bgc.RemOx    = sms.RemOx;
 bgc.Ammox    = sms.Ammox;
 bgc.Anammox  = sms.Anammox;
 bgc.Nitrox    = sms.Nitrox;    
 bgc.RemDen   = sms.RemDen;
 bgc.RemDen1   = sms.RemDen1;
 bgc.RemDen2   = sms.RemDen2;
 bgc.RemDen3   = sms.RemDen3;    
 bgc.Jnn2o_hx = sms.Jnn2o_hx;
 bgc.Jnn2o_nden = sms.Jnn2o_nden;
 bgc.Jno2_hx = sms.Jno2_hx;
 bgc.Jno2_nden = sms.Jno2_nden;
 bgc.Jn2o_prod = sms.Jn2o_prod;
 bgc.Jn2o_cons = sms.Jn2o_cons;
 bgc.Jno2_prod = sms.Jno2_prod;
 bgc.Jno2_cons = sms.Jno2_cons;    
 bgc.sms_n2o = sms.n2o;    

 if bgc.RunIsotopes
 	bgc.r15no3 = sms.r15no3;
 	bgc.r15no2 = sms.r15no2;
 	bgc.r15nh4 = sms.r15nh4;
 	bgc.r15norg = sms.r15nh4;
 	bgc.r15n2o = sms.r15n2o;
 end
