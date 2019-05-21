function bgc = bgc1d_initboundary(bgc)

 if strcmp(bgc.region, 'ETNP')
        bgc.poc_flux_top = -12/86400;   % Bry cond. for top POC Flux            (mmolC/m2/s)
        bgc.o2_top  = 190.0;            % Bry cond. for surface Oxygen          (mmolO2/m3)
        bgc.o2_bot  = 29.;              % Bry cond. for deep Oxygen             (mmolO2/m3)
        bgc.no3_top = 2.9187;           % Bry cond. for surface Nitrate         (mmolNO3/m3)
        bgc.no3_bot = 45.2;             % Bry cond. for deep nitrate            (mmolNO3/m3)
        bgc.po4_top = 0.7450;           % Bry cond. for surface phosphate       (mmolPO4/m3)
        bgc.po4_bot =  3.36;            % Bry cond. for deep phosphate          (mmolPO4/m3)
        bgc.n2o_top = 40/1000;          % Bry cond. for surface nitrous oxide   (mmolN2O/m3)
        bgc.n2o_bot = 39.2/1000;        % Bry cond. for deep nitrous oxide      (mmolN2O/m3)
        bgc.nh4_top = 0.0979;           % Bry cond. for surface ammonia         (mmolNH4/m3)
        bgc.nh4_bot = 0.0111 ;          % Bry cond. for deep ammonia            (mmolNH4/m3)
        bgc.no2_top = 1.0;              % Bry cond. for surface nitrite         (mmolNO2/m3)
        bgc.no2_bot = 0.0411;           % Bry cond. for deep nitrite            (mmolNO2/m3)
        bgc.n2_top = 2.0;               % Bry cond. for surface N2 excess       (mmolN2/m3)
        bgc.n2_bot = 0.8131;            % Bry cond. for deep N2 excess          (mmolN2/m3)
        bgc.d15no3_top = 15.0;          % Bry cond. for surface del15NO3        (permil)
        bgc.d15no3_bot = 5.0;           % Bry cond. for deep del15NO3           (permil)
        bgc.d15no2_top = -10.0;         % Bry cond. for surface del15NO2        (permil)
        bgc.d15no2_bot = -30.0;         % Bry cond. for deep del15NO2           (permil)
        bgc.d15nh4_top = 15.0;          % Bry cond. for surface del15NH4        (permil)
        bgc.d15nh4_bot = 5.0;           % Bry cond. for deep del15NH4           (permil)
        bgc.d15n2oA_top = 7.0;          % Bry cond. for surface delN2O-A        (permil)
        bgc.d15n2oA_bot = 20.0;         % Bry cond. for deep delN2O-A           (permil)
        bgc.d15n2oB_top = -5.0;         % Bry cond. for surface delN2O-B        (permil)
        bgc.d15n2oB_bot = -10.0;        % Bry cond. for deep delN2O-B           (permil)	



  elseif strcmp(bgc.region, 'ETSP')
        bgc.poc_flux_top = -7.0/86400;%-5.96/86400; % Bry cond. for top POC Flux            (mmolC/m2/s)
        bgc.o2_top  = 215;              % Bry cond. for surface Oxygen          (mmolO2/m3)
        bgc.o2_bot  = 77.; %15;%8.2     % Bry cond. for deep Oxygen             (mmolO2/m3)
        bgc.no3_top = 2.813;            % Bry cond. for surface Nitrate         (mmolNO3/m3)
        bgc.no3_bot = 42.5;             % Bry cond. for deep nitrate            (mmolNO3/m3)
        bgc.po4_top = 0.82;             % Bry cond. for surface phosphate       (mmolPO4/m3)
        bgc.po4_bot =  3.055;           % Bry cond. for deep phosphate          (mmolPO4/m3)
        bgc.n2o_top = 13.0/1000;        % Bry cond. for surface nitrous oxide   (mmolN2O/m3)
        bgc.n2o_bot = 35.0/1000;        % Bry cond. for deep nitrous oxide      (mmolN2O/m3)
        bgc.nh4_top = 0.4;              % Bry cond. for surface ammonia         (mmolNH4/m3)
        bgc.nh4_bot = 10^-23 ;             % Bry cond. for deep ammonia            (mmolNH4/m3)
        bgc.no2_top = 0.15;             % Bry cond. for surface nitrite         (mmolNO2/m3)
        bgc.no2_bot = 10^-23;                % Bry cond. for surface nitrite         (mmolNO2/m3)
        bgc.n2_top = 2.0;               % Bry cond. for surface N2 excess       (mmolN2/m3)
        bgc.n2_bot = 6.0;               % Bry cond. for deep N2 excess          (mmolN2/m3)
        bgc.d15no3_top = 15.0;          % Bry cond. for surface del15NO3        (permil)
        bgc.d15no3_bot = 5.0;           % Bry cond. for deep del15NO3           (permil)
        bgc.d15no2_top = 5.0;         % Bry cond. for surface del15NO2        (permil)
        bgc.d15no2_bot = 5.0;         % Bry cond. for deep del15NO2           (permil)
	bgc.d15nh4_top = 7.0;          % Bry cond. for surface del15NH4        (permil)
        bgc.d15nh4_bot = 7.0;           % Bry cond. for deep del15NH4           (permil)
        bgc.d15n2oA_top = 7.0;          % Bry cond. for surface delN2O-A        (permil)
        bgc.d15n2oA_bot = 20.0;         % Bry cond. for deep delN2O-A           (permil)
        bgc.d15n2oB_top = -5.0;         % Bry cond. for surface delN2O-B        (permil)
        bgc.d15n2oB_bot = -10.0;        % Bry cond. for deep delN2O-B           (permil)

  end
