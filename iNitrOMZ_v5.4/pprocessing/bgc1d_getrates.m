function bgc = bgc1d_getrates(bgc,Data)
	zstep = bgc.npt;
	bgc.z = linspace(bgc.zbottom,bgc.ztop,zstep+1);
	bgc.dz = (bgc.ztop - bgc.zbottom) / zstep;
	bgc.SolNames = {'o2','no3','poc','po4','n2o', 'nh4', 'no2', 'n2'};
	ntrSol = length(bgc.SolNames);
	for indt=1:ntrSol
		bgc.(bgc.SolNames{indt}) = bgc.sol(indt,:);
	end
	sms =  bgc1d_sms_diag(bgc);
	bgc.rates = nan(length(Data.rates.name), length(bgc.zgrid));
        tmp.nh4ton2o = sms.SoN2OAo;
        tmp.noxton2o = bgc.NCden2*sms.RemDen2;
        tmp.nh4tono2 = sms.Ammox;
        tmp.no3tono2 = bgc.NCden1*sms.RemDen1;
        tmp.anammox = sms.Anammox;
	for n = 1 : length(Data.rates.name)
		bgc.rates(n,:) = tmp.(Data.rates.name{n});
	end
