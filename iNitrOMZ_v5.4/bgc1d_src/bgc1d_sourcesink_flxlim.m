 function sms =  bgc1d_sourcesink(bgc,t); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % iNitrOMZ v1.0 - Simon Yang  - April 2019
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Specifies the biogeochemical sources and sinks for OMZ nutrient/POC/N2O model
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Make sure we have no negative concentration
% ff=fields(t);
% for f = 1 : length(ff)
% 	t.(ff{f})(t.(ff{f})<0)=0;
% end

 % % % % % % % % % % % %
 % % % % J-OXIC  % % % %
 % % % % % % % % % % % %

 %----------------------------------------------------------------------
 % (1) Oxic Respiration rate (C-units):
 %----------------------------------------------------------------------
 RemOx = bgc.Krem .* mm1(t.o2,bgc.KO2Rem) .* t.poc;

 %----------------------------------------------------------------------
 % (2) Ammonium oxidation (molN-units):
 %----------------------------------------------------------------------
 Ammox = bgc.KAo .*  mm1(t.o2,bgc.KO2Ao) .*  mm1(t.nh4,bgc.KNH4Ao) ;

 %----------------------------------------------------------------------
 % (3) Nitrite oxidation (molN-units):
 %----------------------------------------------------------------------
 Nitrox = bgc.KNo .*  mm1(t.o2,bgc.KO2No) .* mm1(t.no2,bgc.KNO2No);

 %----------------------------------------------------------------------
 % (4) N2O and NO2 production by ammox and nitrifier-denitrif (molN-units): 
 %----------------------------------------------------------------------
 Y = n2o_yield(t.o2, bgc);
 % via NH2OH
 Jnn2o_hx = Ammox .* Y.nn2o_hx_nh4;
 Jno2_hx = Ammox .* Y.no2_hx_nh4;;
 % via NH4->NO2->N2O
 Jnn2o_nden = Ammox .* Y.nn2o_nden_nh4;
 Jno2_nden = Ammox .* Y.no2_nden_nh4;
 
 % % % % % % % % % % % %
 % % %   J-ANOXIC  % % % 
 % % % % % % % % % % % %

 %----------------------------------------------------------------------
 % (5) Denitrification (C-units)
 %----------------------------------------------------------------------
 RemDen1 = bgc.KDen1.* mm1(t.no3,bgc.KNO3Den1) .* fexp(t.o2,bgc.KO2Den1) .* t.poc;
 RemDen2 = bgc.KDen2 .* mm1(t.no2,bgc.KNO2Den2) .* fexp(t.o2,bgc.KO2Den2) .* t.poc;
 RemDen3 = bgc.KDen3 .* mm1(t.n2o,bgc.KN2ODen3) .* fexp(t.o2,bgc.KO2Den3) .* t.poc;

 %----------------------------------------------------------------------
 % (6) Anaerobic ammonium oxidation (molN-units):
 %----------------------------------------------------------------------
 Anammox = bgc.KAx .* mm1(t.nh4,bgc.KNH4Ax) .* mm1(t.no2,bgc.KNO2Ax) .* fexp(t.o2,bgc.KO2Ax);

 %----------------------------------------------------------------------
 % (7)  Prevent overshooting to negative concentrations
 %----------------------------------------------------------------------
 if (t.o2 +(-bgc.OCrem .* RemOx - 1.5*Ammox - 0.5 .* Nitrox).*bgc.dt < 0)
 	RemOx = (t.o2./bgc.dt -1.5*Ammox - 0.5 .* Nitrox)./bgc.OCrem;
 	RemOx(RemOx<0)=0;
 elseif t.no3+(-bgc.NCden1 .* RemDen1 + Nitrox).*bgc.dt < 0
 	RemDen1 = (t.no3./bgc.dt + Nitrox)./bgc.NCden1;	
 	RemDen1(RemDen1<0)=0;
 elseif t.n2o+( 0.5 * bgc.NCden2 * RemDen2 - bgc.NCden3 * RemDen3 + 0.5 .* (Jnn2o_hx+Jnn2o_nden)).*bgc.dt<0
 	RemDen3 = (t.n2o./bgc.dt + 0.5 * bgc.NCden2 * RemDen2 +  0.5 .* (Jnn2o_hx+Jnn2o_nden))./ bgc.NCden3;
 	RemDen3(RemDen3<0)=0;
 elseif t.no2+(Jno2_hx + Jno2_nden - Nitrox + bgc.NCden1 .* RemDen1 - bgc.NCden2 .* RemDen2 - Anammox).*bgc.dt<0
 	RemDen2 = (+t.no2./bgc.dt+Jno2_hx + Jno2_nden - Nitrox + bgc.NCden1 .* RemDen1 - Anammox)./bgc.NCden2;
 	RemDen2(RemDen2<0)=0;
 end
% % rescale remin. rates 
 KRemOx = RemOx./t.poc;
 KRemDen1 = RemDen1./t.poc;
 KRemDen2 = RemDen2./t.poc;
 KRemDen3 = RemDen3./t.poc;
 
 %----------------------------------------------------------------------
 % (8)  Calculate SMS for each tracer
 %---------------------------------------------------------------------- 
 sms.o2   =  (-bgc.OCrem .* RemOx - 1.5.*Ammox - 0.5 .* Nitrox);
 sms.no3  =  (Nitrox - bgc.NCden1 .* RemDen1);
 sms.poc  =  (- RemOx - (RemDen1+RemDen2+RemDen3));
 sms.po4  =  (+bgc.PCrem .* (RemOx) + bgc.PCden1 .* RemDen1 + bgc.PCden2 .* RemDen2 + bgc.PCden3 .* RemDen3);
 sms.nh4  =  (+ bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3) - Ammox - Anammox);
 sms.no2  =  (Jno2_hx + Jno2_nden + bgc.NCden1 .* RemDen1 - bgc.NCden2 .* RemDen2 - Anammox - Nitrox);
 sms.n2   =  (bgc.NCden3 .* RemDen3 + Anammox);
 sms.kpoc = -(KRemOx + KRemDen1 + KRemDen2 + KRemDen3);
 % N2O individual SMSs
 sms.n2oind.ammox = 0.5 .* Jnn2o_hx;
 sms.n2oind.nden  = 0.5 .* Jnn2o_nden;
 sms.n2oind.den2  = 0.5 .* bgc.NCden2 .* RemDen2;
 sms.n2oind.den3  = - bgc.NCden3 .* RemDen3;
 % N2O total SMS
 sms.n2o = sms.n2oind.ammox + sms.n2oind.nden + sms.n2oind.den2 + sms.n2oind.den3;

 if bgc.RunIsotopes
	 % update 15N/N ratios
	 bgc = bgc1d_initIso_update_r15n(bgc,t);
	 % Calculate sources and sinks for 15N tracers
	 sms.i15no3 = bgc.r15no2 .* bgc.alpha_nitrox .* Nitrox ...
	 	    - bgc.r15no3 .* bgc.alpha_den1 .* bgc.NCden1 .* RemDen1;
	 sms.i15no2 = bgc.r15nh4 .* (bgc.alpha_ammox_no2 .* Jno2_hx + bgc.alpha_nden_no2 .* Jno2_nden) ...
	            + bgc.r15no3 .* bgc.alpha_den1 .* bgc.NCden1 .* RemDen1 ...
	    	    - bgc.r15no2 .* bgc.alpha_den2 .* bgc.NCden2 .* RemDen2 ...
	            - bgc.r15no2 .* bgc.alpha_ax_no2 .* Anammox ...
	            - bgc.r15no2 .* bgc.alpha_nitrox .* Nitrox ;
	 sms.i15nh4 = bgc.r15norg.* bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3) ...
	            - bgc.r15nh4 .* (bgc.alpha_ammox_no2 .* Jno2_hx + bgc.alpha_nden_no2 .* Jno2_nden) ...
	 	    - bgc.r15nh4 .* (bgc.alpha_ammox_n2o .* Jnn2o_hx + bgc.alpha_nden_n2o .* Jnn2o_nden) ...
		    - bgc.r15nh4 .* bgc.alpha_ax_nh4 .* Anammox;
	 % N2O indivisual SMS	    
	 sms.i15n2oind.ammox =  bgc.r15nh4 .* bgc.alpha_ammox_n2o .* sms.n2oind.ammox;
	 sms.i15n2oind.nden  =  bgc.r15nh4 .* bgc.alpha_nden_n2o .* sms.n2oind.nden;
  	 sms.i15n2oind.den2  =  bgc.r15no2 .* bgc.alpha_den2 .* sms.n2oind.den2;
	 % Get isotopomer partitionning
	 %ammox
	 ii = dNisoSP('i15N', sms.i15n2oind.ammox, 'N', sms.n2oind.ammox, 'SP', bgc.n2oSP_ammox); 
	 sms.i15n2oAind.ammox = ii.i15N_A;
	 sms.i15n2oBind.ammox = ii.i15N_B;
	 %nitrifier-denitrification
	 ii = dNisoSP('i15N', sms.i15n2oind.nden, 'N', sms.n2oind.nden, 'SP', bgc.n2oSP_nden);
	 sms.i15n2oAind.nden = ii.i15N_A;
         sms.i15n2oBind.nden = ii.i15N_B;
	 %denitrification 2 (no2-->n2o)
	 ii = dNisoSP('i15N', sms.i15n2oind.den2, 'N', sms.n2oind.den2, 'SP', bgc.n2oSP_den2);
	 sms.i15n2oAind.den2 = ii.i15N_A;
         sms.i15n2oBind.den2 = ii.i15N_B;
	 % Total
	 sms.i15n2oA = ii.i15N_A + ii.i15N_A + ii.i15N_A + bgc.r15n2oA .* bgc.alpha_den3_Alpha .* sms.n2oind.den3;
	 sms.i15n2oB = ii.i15N_B + ii.i15N_B + ii.i15N_B + bgc.r15n2oB .* bgc.alpha_den3_Beta .* sms.n2oind.den3;
	 sms.i15n2oind.den3  = bgc.r15n2oA .* bgc.alpha_den3_Alpha .* sms.n2oind.den3 + bgc.r15n2oB .* bgc.alpha_den3_Beta .* sms.n2oind.den3;
         % Sum all SMS
         sms.i15n2o = sms.i15n2oind.ammox + sms.i15n2oind.nden + sms.i15n2oind.den2 + sms.i15n2oind.den3;

	 % Checks
	 if (sms.i15n2oind.ammox ~= sms.i15n2oAind.ammox + sms.i15n2oBind.ammox) 
	    error('Ammox: prod. of i15n2o is not equal to the summed prod. of i15n2oA and i15n2oB')
    	 end
	 if (sms.i15n2oind.nden ~= sms.i15n2oAind.nden + sms.i15n2oBind.nden)
            error('Nitrifier-denitrificaiton: prod. of i15n2o is not equal to the summed prod. of i15n2oA and i15n2oB')
         end
	 if (sms.i15n2oind.den2 ~= sms.i15n2oAind.den2 + sms.i15n2oBind.den2)
            error('Nitrite reduction (den2): prod. of i15n2o is not equal to the summed prod. of i15n2oA and i15n2oB')
         end

 end

