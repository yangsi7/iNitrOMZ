function bgc = bgc1d_initIso_params(bgc)


 %%%%%%%% Ammonification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgc.d15norg = 7.0;               % d15N of organic matter (permil)


 %%%%%% Ammonium oxidation %%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Epsilon for ...
 % .. ammox -- NH4-->NO2: 
 bgc.eps_ammox_no2 = 15.0;          % Enrichment factors during ammox (permil)
 % .. ammox -- NH4-->N2O:
 bgc.eps_ammox_n2o = 47.0;          % Enrichment factors during N2O prod. via NH2OH (permil)
 % .. nitrif-denitrif -- NH4-->NO2:
 bgc.eps_nden_no2 = 15.0;           % Enrichment factors during nitrifier-denitrification (permil)
 % .. nitrif-denitrif -- NH4-->N2O:
 bgc.eps_nden_n2o = 58.0;           % Enrichment factors during nitrifier-denitrification (permil) 

 % SP for ...
 % .. ammox -- NH4-->N2O:
 bgc.n2oSP_ammox = 36.0;           % Site preference during N2O prod. via NH2OH (permil)
 % .. nitrif-denitrif -- NH4-->N2O:
 bgc.n2oSP_nden = -11.0;           % Site preference during nitrifier-denitrification (permil)


 %%%%%% Nitrite oxidation %%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Epsilon for nitrox -- NO2-->NO3
 bgc.eps_nitrox = -12.8;          % Enrichment factors during nitrox (permil)   

 %%%%%%% Denitrification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % epsilon for den1 -- NO3-->NO2
 bgc.eps_den1 = 30.0;          % Enrichment factors during denitrification 1 (permil)
 % epsilon for den2 -- NO2-->N2O
 bgc.eps_den2 = 25.0;          % Enrichment factors during denitrification 2 (permil)
 % epsilon for den3 -- N2O-->N2
 bgc.eps_den3 = 13.0;          % Enrichment factors during denitrification 3 (permil)

 % SP for ..
 % .. den2 -- NO2-->N2O:
 bgc.n2oSP_den2 = -1.0;           % Site preference during denitrification 2 (permil)

 %%%%%%%%%%% Anammox %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Epsilon for Anammox -- NO2-->N2
 bgc.eps_ax_no2 = -40.0;          % Enrichment factors during anammox (permil)
 % Epsilon for Anammox -- NH4-->N2
 bgc.eps_ax_nh4 = 25.0;          % Enrichment factors during anammox (permil)

 bgc = bgc1d_initIso_Dep_params(bgc);
