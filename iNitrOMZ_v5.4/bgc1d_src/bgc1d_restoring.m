 function restoring = bgc1d_restoring(bgc,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d model v1.3 - Simon Yang  - October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies restoring source (linear relaxation) for selected dissolved tracers
 if bgc.tauZvar == 1
    tauh = ppval(bgc.restore.tauh_pre,bgc.zgrid);
 else
    tauh = bgc.tauh;
 end
 
 if bgc.PO4rest==0
    restoring.po4 = zeros(1,length(t.po4));
 else
    restoring.po4 = (bgc.restore.po4_cout-t.po4)./tauh;
 end

 if bgc.NO3rest==0
    restoring.no3 = zeros(1,length(t.no3));
 else
    restoring.no3 = (bgc.restore.no3_cout-t.no3)./tauh; 
 end

 if bgc.O2rest==0
    restoring.o2 = zeros(1,length(t.o2));
 else
    if bgc.forceanoxic == 1
	cout(find(bgc.zgrid==bgc.forceanoxic_bounds(2)):find(bgc.zgrid==bgc.forceanoxic_bounds(1)))=0;
    end
    restoring.o2 = (bgc.restore.o2_cout-t.o2)./tauh;
 end

 if bgc.N2Orest==0
    restoring.n2o = zeros(1,length(t.n2o));
 else
    restoring.n2o = (bgc.restore.n2o_cout-t.n2o)./tauh;
 end

 if bgc.NO2rest==0
    restoring.no2 = zeros(1,length(t.no2));
 else
    restoring.no2 = (bgc.restore.no2_cout-t.no2)./tauh;
 end

 if bgc.NH4rest==0
    restoring.nh4 = zeros(1,length(t.nh4));
 else
    restoring.nh4 = (bgc.restore.nh4_cout-t.nh4)./tauh;
 end

 if bgc.N2rest==0
    restoring.n2 = zeros(1,length(t.n2));
 else
    cout = zeros(size(bgc.zgrid));
    restoring.n2 = (bgc.restore.n2_cout-t.n2)./tauh;
 end
 if bgc.RunIsotopes
 	if bgc.i15NO3rest==0
 	   restoring.i15no3 = zeros(1,length(t.i15no3));
 	else
 	   cout = zeros(size(bgc.zgrid));
 	   restoring.i15no3 = (bgc.restore.i15no3_cout-t.i15no3)./tauh;
 	end

 	if bgc.i15NO2rest==0
 	   restoring.i15no2 = zeros(1,length(t.i15no2));
 	else
 	   cout = zeros(size(bgc.zgrid));
 	   restoring.i15no2 = (bgc.restore.i15no2_cout-t.i15no2)./tauh;
 	end

 	if bgc.i15NH4rest==0
 	   restoring.i15nh4 = zeros(1,length(t.i15nh4));
 	else
 	   cout = zeros(size(bgc.zgrid));
 	   restoring.i15nh4 = (bgc.restore.i15nh4_cout-t.i15nh4)./tauh;
 	end 

 	if bgc.i15N2OArest==0
 	   restoring.i15n2oA = zeros(1,length(t.i15n2oA));
 	else
 	   cout = zeros(size(bgc.zgrid));
 	   restoring.i15n2oA = (bgc.restore.i15no2A_cout-t.i15n2oA)./tauh;
 	end

 	if bgc.i15N2OBrest==0
 	   restoring.i15n2oB = zeros(1,length(t.i15n2oB));
 	else
 	   cout = zeros(size(bgc.zgrid));
 	   restoring.i15n2oB = (bgc.restore.i15n2ob_cout-t.i15n2oB)./tauh;
 	end 
 end

