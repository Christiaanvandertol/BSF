% first run the algorithm. Then we have a lot in the workshpace to work
% with.

%master_selecteddays_revision;
%master_Cabauw2;

%% the FWMH values for which the effect of the spectral response function is evaluated
FWHMi       = (.01:.01:.7)';

%% parameters for the retrieval algorithm
stoptol     = 1E-6;
opt         = optimset('MaxIter',30,'TolFun',stoptol);
aprior      = 1;

%% loading SCOPE reflectances and fluorescence
SCOPEoutputfolder = 'OHP_2021-06-25-1444';
rsd         = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\rsd.csv'],',',2,0); %#ok<*DLMRD>
rdd         = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\rdd.csv'],',',2,0);
fscope      = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\fluorescence.csv'],',',2,0);
wlS         = load(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\wlS.txt']);

%% the wavelength of the FLoX (modify this statement for your system)
wl          = D(1).wl;

%%
SZA         = [75,60,45,30];
cos_sza     = cos(SZA/180*pi);

%% load MODTRAN data
%soltir_tp7([ModtranFile{site} '.tp7']);
s1          = importdata('FLEX-S3_SPAIN_SZA75.atm');
s2          = importdata('FLEX-S3_SPAIN_SZA60.atm');
s3          = importdata('FLEX-S3_SPAIN_SZA45.atm');
s4          = importdata('FLEX-S3_SPAIN_SZA30.atm');

%% Experiment 1. For a fixed value of FWMH (of 0.31), error as function of atmospheric path length
% In this piece of code, we also estimate a correction function for the
% effect of the SRF on the spectral shape of the O2A and O2B bands. This
% is estimated as follows:
% We estimate <L> as:
% Second: we estimate L_{approx} as: exp(log(<E>r)*a
% Third: SRC = (<L>-L_{approx}) / (1-a)

ai          = (1.001:.004:1.1);
priorweight = 0;
FWHM        = 0.31;
sigma       = FWHM/2.355; % see https://en.wikipedia.org/wiki/Full_width_at_half_maximum

[piL,E]     = deal(NaN*wl);
[ak,Fk,FFLDk,EX]= deal(NaN*ones(length(ai),1));
f           = fscope(5,760-639);
nwl         = [59,9];

for O2band = 1:2
    SRC_SZA     = zeros(nwl(O2band),4);
    SRCj        = zeros(nwl(O2band),25,4);

    for iSZA = 1:length(SZA)

        switch iSZA
            case 1, s = s1;
            case 2, s = s2;
            case 3, s = s3;
            case 4, s = s4;
        end

        % interpolations and pre-processing
        wl_MODTRAN  = transform_wvl_from_vac_to_air(s.data(:,2));


        T           = s.data(:,3:20);
        atmo.M      = [T(:,1) T(:,3) T(:,4) T(:,5) T(:,12) T(:,16)];
        atmo.Ta     = 25; % influence of air temperature is negligible, but a value
        % has to be provided.

        %interpolate SCOPE simulated reflectance to MODTRAN wavelengths
        SAIL.rsd    = interp1(wlS,rsd(5,:),wl_MODTRAN);
        SAIL.rdd    = interp1(wlS,rdd(5,:),wl_MODTRAN);

        % calculate irradiance
        [Esun,Esky] = calcIrradiance(atmo,SAIL,wl_MODTRAN);
        % Esky(isnan(Esky))=0;

        % interpolate SCOPE fluorescence to the MODTRAN wavelengths
        F_MODTRAN   = interp1((640:850),fscope(5,:),wl_MODTRAN) * sum(Esun)/8.3980e+06; %scale fluorescence
        E_MODTRAN   = Esun;%+Esky;
        SAIL.rsd(isnan(SAIL.rsd)) = 0;
        F_MODTRAN(isnan(F_MODTRAN))= 0;

        % normalize the band depth by the interpolated values
        piL_MODTRAN0    = E_MODTRAN.*SAIL.rsd;

        index2      = find(wl_MODTRAN>p.wl_left(O2band) & wl_MODTRAN<p.wl_right(O2band));
        normpiL2    = interp1([wl_MODTRAN(index2(1)),wl_MODTRAN(index2(end))],[piL_MODTRAN0(index2(1)),piL_MODTRAN0(index2(end))],wl_MODTRAN(index2));

        for j = 1:length(ai)
            a = ai(j);
            aprior = a;

            piL_MODTRAN         = piL_MODTRAN0;
            piL_MODTRAN(index2) = normpiL2.* exp(log(piL_MODTRAN(index2)./normpiL2)*a);

            % convolution
            for k = 1:length(wl)
                y = normpdf(wl_MODTRAN,wl(k),sigma);
                E(k) = sum(y.*E_MODTRAN)/sum(y);
                piL(k) = sum(y.*piL_MODTRAN)/sum(y);
            end

            [O2A,O2B] = retrievalF(wl,E,piL,opt,a,0.7,1,priorweight,p,zeros(nwl(1),1),zeros(nwl(2),1));

            switch O2band
                case 1, SRCj(:,j,iSZA) = (log(O2A.piL./O2A.normpiL) - aprior* log(O2A.E./O2A.normE));
                case 2, SRCj(:,j,iSZA) = (log(O2B.piL./O2B.normpiL) - aprior* log(O2B.E./O2B.normE));
            end

            for k = 1:size(SRCj,1)
                SRC_SZA(k,iSZA) = (ai-1)' \ SRCj(k,:,iSZA)';
            end

  %          if j == 7 && iSZA==3
  %               switch O2band
  %                  case 1, O2 = O2A;
  %                  case 2, O2 = O2B;
  %              end
%                 figure(2)
%                 subplot(2,2,O2band*2-1)
%                 plot(O2.wl,[(log(O2.piL./O2.normpiL)), aprior* log(O2.E./O2.normE)]), hold on
%                 xlabel('wl (nm)')
%                 ylabel('log(\piL/\piL_{norm})')
% 
%                 subplot(2,2,2*O2band)
%                 plot(O2.wl,(log(O2.piL./O2.normpiL))-aprior* log(O2.E./O2.normE))
%                 xlabel('wl (nm)')
%                 ylabel('SRC')
% 
%             end


        end

    end

    SRC = SRC_SZA(:,4)./mean(SRC_SZA(:,4));
    c = polyfit((SZA), mean(SRCA_SZA),3);

    switch O2band
        case 1
            SRCA = SRC;
            %save('SRCA_ground.mat','SRCA', 'c')
            save('SRCA.mat','SRCA', 'c')
        case 2
            SRCB = SRC;
            %save('SRCB_ground.mat','SRCB', 'c')
            save('SRCB.mat','SRCB', 'c')
    end

%     figure(1)
%     subplot(2,1,O2band)
%     plot(SZA,mean(SRC_SZA),'kx')
%     hold on
%     plot([30:75],polyval(c,[30:75]),'k')
%     xlabel('SZA')
%     ylabel('mean SRC/(1-a)')

% figure(3), hold on
% meanSRCj = mean(SRCj,[1,3]);
% plot(ai,meanSRCj,'x')
end


%%
figure(3), hold on
meanSRCj = mean(SRCj,[1,3]);
plot(ai,meanSRCj,'x')
xlabel('a')
ylabel('SRC_{mean}')
