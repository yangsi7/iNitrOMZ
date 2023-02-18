function bgc1d_plotisotopes(bgc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d ncycle v 1.0 - Simon Yang  - October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results from optimization by the genetic algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot modelled O2, NO3, PO4, POC, N2O vs the Data used for the optimization%%%

 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 10 7])

 idx = find(bgc.no2>0.005);
 blue = [0, 0.4470, 0.7410];
 red = [0.8500, 0.3250, 0.0980];
 yellow = [0.9290, 0.6940, 0.1250];
 purple = [0.4940, 0.1840, 0.5560];
 green = [0.4660, 0.6740, 0.1880];

 subplot(2,2,1)
 plot(bgc.d15no3,bgc.zgrid,'color', blue,'linewidth',3, 'DisplayName','\delta^{15}NO_3^-')
 grid on; hold on
 plot(bgc.d15no2(idx),bgc.zgrid(idx),'color', red,'linewidth',3, 'DisplayName','d15NO2-')
 grid on; hold on
 plot(bgc.d15nh4,bgc.zgrid,'color', yellow, 'linewidth',3, 'DisplayName','d15NH4+')
 title('Substrates')
 ylabel('z (m)')
 xlabel(char(8240))
 ylim( [bgc.zbottom bgc.ztop]);
 %xlim( [-5 max(bgc.o2)]);
 legend('\delta^{15}N-NO_3^-', '\delta^{15}N-NO_2^-', '\delta^{15}N-NH_4^+', 'Location', 'southwest')
 grid off; hold off

 subplot(2,2,2)
 plot(bgc.n2o*1000,bgc.zgrid,'k','linewidth',3)
 title('[N2O]')
 ylabel('z (m)')
 xlabel('[N_2O] (nM)')
 ylim( [bgc.zbottom bgc.ztop]);
 %xlim( [-5 max(bgc.o2)]);

 subplot(2,2,3)
 plot(bgc.d15n2oA,bgc.zgrid,'color', purple,'linewidth',3, 'DisplayName','d15n2oA')
 grid on; hold on
 plot(bgc.d15n2oB,bgc.zgrid,'color', green,'linewidth',3, 'DisplayName','d15n2oB')
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = [0 0.5 0.5];
 title('d15N2O')
 ylabel('z (m)')
 xlabel(char(8240))
 legend('\delta^{15}N-N_2O^{\alpha}','\delta^{15}N-N_2O^{\beta}','Location','southeast')
 ylim( [bgc.zbottom bgc.ztop]);

 subplot(2,2,4)
 plot((bgc.d15n2oA - bgc.d15n2oB),bgc.zgrid,'k','linewidth',3)
 title('SP')
 ylabel('z (m)')
 xlabel(char(8240))
 ylim( [bgc.zbottom bgc.ztop]);
 %xlim( [-5 max(bgc.o2)]);

 end

