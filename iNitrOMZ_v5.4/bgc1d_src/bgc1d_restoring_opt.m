 function restoring = bgc1d_restoring(bgc,Soltime);
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
    restoring(4,:) = 0;
 else
    cout = ppval(bgc.restore.po4_cs,bgc.zgrid);
    restoring(4,:) = (cout-Soltime(4,:))./tauh;
 end
 if bgc.NO3rest==0
    restoring(3,:) = 0;
 else
    cout = ppval(bgc.restore.no3_cs,bgc.zgrid);
    restoring(3,:) = (cout-Soltime(3,:))./tauh; 
 end
 if bgc.O2rest==0
    restoring(2,:) = 0;
 else
    cout = ppval(bgc.restore.o2_cs,bgc.zgrid);
    if bgc.forceanoxic == 1
	cout(find(bgc.zgrid==bgc.forceanoxic_bounds(2)):find(bgc.zgrid==bgc.forceanoxic_bounds(1)))=0;
    end
    restoring(2,:) = (cout-Soltime(2,:))./tauh;
 end
 if bgc.N2Orest==0
    restoring(5,:) = 0;
 else
    cout = ppval(bgc.restore.n2o_cs,bgc.zgrid);
    restoring(5,:) = (cout-Soltime(5,:))./tauh;
 end
  if bgc.NO2rest==0
    restoring(7,:) = 0;
 else
    cout = ppval(bgc.restore.no2_cs,bgc.zgrid);
    restoring(7,:) = (cout-Soltime(7,:))./tauh;
  end
   if bgc.NH4rest==0
    restoring(6,:) = 0;
 else
    cout = ppval(bgc.restore.nh4_cs,bgc.zgrid);
    restoring(6,:) = (cout-Soltime(6,:))./tauh;
 end
    if bgc.N2rest==0
    restoring(8,:) = 0;
 else
    cout = zeros(size(bgc.zgrid));
    restoring(8,:) = (cout-Soltime(8,:))./tauh;
 end

