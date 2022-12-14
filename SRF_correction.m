% first run the algorithm. Then we have a lot in the workshpace to work
% with.

%master_selecteddays_revision;

%% options
experiment1 = 1;
experiment2 = 0;
experiment3 = 0;

%% parameters for plotting of the results
alphabet = {'a','b','c','d'};
colors      = [0 0 1;1 0 0];

%% the FWMH values for which the effect of the spectral response function is evaluated
FWHMi       = (.01:.01:.7)';

%% parameters for the retrieval algorithm
stoptol     = 1E-6;
opt         = optimset('MaxIter',30,'TolFun',stoptol);

%% loading SCOPE reflectances and fluorescence
SCOPEoutputfolder = 'OHP_2021-06-25-1444';
rsd         = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\rsd.csv'],',',2,0); %#ok<*DLMRD>
rdd         = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\rdd.csv'],',',2,0);
fscope      = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\fluorescence.csv'],',',2,0);
wlS         = load(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\wlS.txt']);

%% loading (earlier derived, with this script) correction spectra
%load SRCA.mat
%load SRCB.mat

%% loading (earlier derived, with this script) correction spectra
SRCA_SZA    = zeros(59,4);
SRCB_SZA    = zeros(10,4);
SRCAj       = zeros(59,25,4);
SRCBj       = zeros(10,25,4);

%% the wavelength of the FLoX
wl          = D(1).wl;

%% load MODTRAN data
%soltir_tp7([ModtranFile{site} '.tp7']);
s1          = importdata('FLEX-S3_SPAIN_SZA75.atm');
s2          = importdata('FLEX-S3_SPAIN_SZA60.atm');
s3          = importdata('FLEX-S3_SPAIN_SZA45.atm');
s4          = importdata('FLEX-S3_SPAIN_SZA30.atm');
%s1          = importdata('FLEX-S3_SPAIN_SZA30.atm');
%s2          = importdata('FLEX-S3_GER1_PC.atm');

%% Experiment 1. For a fixed value of FWMH (of 0.31), error as function of atmospheric path length
% In this piece of code, we also estimate a correction function for the
% effect of the SRF on the spectral shape of the O2A and O2B bands. This
% is estimated as follows:
% We estimate <L> as:
% Second: we estimate L_{approx} as: exp(log(<E>r)*a
% Third: SRC = (<L>-L_{approx}) / (1-a)

ai          = (1.004:.004:1.1);
priorweight = 0;
FWHM        = 0.31;
sigma       = FWHM/2.355; % see https://en.wikipedia.org/wiki/Full_width_at_half_maximum

[piL,E]     = deal(NaN*wl);
[ak,Fk,FFLDk,EX]= deal(NaN*ones(length(ai),1));
f           = fscope(5,760-639);

SZA         = [30,45,60,75];
cos_zsa     = cos(SZA/180*pi);

for iSZA = 1:length(SZA)
    switch iSZA
        case 1, s = s1;
        case 2, s = s2;
        case 3, s = s3;
        case 4, s = s4;
    end

    % interpolations and pre-processing
    wl_MODTRAN  = transform_wvl_from_vac_to_air(s.data(:,2));
    index2      = find(wl_MODTRAN>p.wl_left(1) & wl_MODTRAN<p.wl_right(1));

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
    normpiL         = interp1([wl_MODTRAN(index2(1)),wl_MODTRAN(index2(end))],[piL_MODTRAN0(index2(1)),piL_MODTRAN0(index2(end))],wl_MODTRAN(index2));

    for b = 1:3
        % three cases: b=1: no SIF, no correction; b=2: SIF, no correction;
        % b=3: SIF+ correction
        for j = 1:length(ai)
            a = ai(j);
            aprior = a;

            piL_MODTRAN     = piL_MODTRAN0;
            piL_MODTRAN(index2) = normpiL.* exp(log(piL_MODTRAN(index2)./normpiL)*a);

            if b > 1
                piL_MODTRAN = piL_MODTRAN+ F_MODTRAN;
            end

            % convolution
            for k = 1:length(wl)
                y = normpdf(wl_MODTRAN,wl(k),sigma);
                E(k) = sum(y.*E_MODTRAN)/sum(y);
                piL(k) = sum(y.*piL_MODTRAN)/sum(y);
            end

            [O2A,O2B] = retrievalF(wl,E,piL,opt,aprior,0.7,1,priorweight,p,SRCA/cos_sza(SZA(iSZA))*(b==3),SRCB/cos_sza(SZA(iSZA))*(b==3));

            Fk(j)       = O2A.F;
            FFLDk(j)    = O2A.iFLD;
            ak(j)       = O2A.a;
            EX(j)       = O2A.EXITFLAG;

            if b==1
                SRCAj(:,j,iSZA) = (log(O2A.piL./O2A.normpiL) - aprior* log(O2A.E./O2A.normE));
                SRCBj(:,j,iSZA) = (log(O2B.piL./O2B.normpiL) - aprior* log(O2B.E./O2B.normE));
            end


            if iSZA==4 && j==6 && b<3 % just for one value of atmospheric path length (1.02) and SZA of 30 deg
                figure(10), hold on
                subplot(1,3,1)
                plot(O2A.wl,log(O2A.piL./O2A.normpiL) - aprior*log(O2A.E./O2A.normE),'Color',colors(b,:));
                hold on
                ylabel('residual = log(L_{norm}) - a log(E_{norm})')
                xlabel('wl (nm)')
                set(gca, 'FontSize',12)

                subplot(1,3,2) % for O3B
                plot(O2B.wl,log(O2B.piL./O2B.normpiL) - aprior*log(O2B.E./O2B.normE),'Color',colors(b,:));%*(j+2)/(length(ai)+2))
                %set(gca, 'ylim',[0, 0.03])
                ylabel('residual = log(L_{norm}) - a log(E_{norm})')
                xlabel('wl (nm)')
                set(gca, 'FontSize',12)
                hold on
            end
            if iSZA==2 &&j==6 && b==1 %extra fig: for SZA of 60 deg
                figure(10)
                subplot(1,2,1)
                plot(O2A.wl,log(O2A.piL./O2A.normpiL) - aprior*log(O2A.E./O2A.normE),'Color','c');
                hold on
                subplot(1,2,2)
                plot(O2B.wl,log(O2B.piL./O2B.normpiL) - aprior*log(O2B.E./O2B.normE),'Color','c');
                hold on
            end

        end
        if b == 1
            for k = 1:size(SRCAj,1)
                SRCA_SZA(k,iSZA) = (ai-1)' \ SRCAj(k,:,iSZA)';
            end
            for k = 1:size(SRCBj,1)
                SRCB_SZA(k,iSZA) = (ai-1)' \ SRCBj(k,:,iSZA)';
            end
        end

        if b>1 && iSZA==4 % with SIF, and without and with correction for the SRF
            figure(8)
            subplot(2,1,1), z=plot(ai,ak,'kx');
            if b == 2, set(z, 'MarkerEdgeColor','b'), end
            hold on
            plot([1 1.1],[1,1.1],'k')
            ylabel('a_{retrieved}')

            subplot(2,1,2), z = plot(ai,Fk,'kx'); hold on
            if b == 2, set(z, 'MarkerEdgeColor','b'), end
            plot([1 1.1],[f,f],'k')
            ylabel('SIF (Wm^{-1}\mum^{-1}sr^{-1})')
            xlabel('a_{input}')
        end
    end
end

SRCA = SRCA_SZA(:,4)*cos(30/180*pi);
SRCB = SRCB_SZA(:,4)*cos(30/180*pi);
%%
figure(10), hold on

subplot(131)
plot(O2A.wl,SRCAj(:,6,2),'c')
legend('\theta_s=30^o','\theta_s=30^o, with F','\theta_s=60^o')
title('a')

subplot(132), 
title('b')
subplot(133)
cos_zsa = cos([75,60,45,30]/180*pi);
for k = 1:4, plot(ai,SRCAj(14,:,k)*cos_zsa(k),'x'), hold on, end
ylabel('(log(L_{norm}) - a log(E_{norm}))*cos(\theta_s)')
xlabel('a')
set(gca, 'FontSize',12)
title('c')
legend('\theta_s = 75^o', '\theta_s = 60^o','\theta_s = 45^o','\theta_s = 30^o', 'Location','northwest')

%%
figure(100), 
%subplot(1,2,1)
%   plot(O2A.wl,log(O2A.piL./O2A.normpiL) - aprior*log(O2A.E./O2A.normE),'Color',colors(b,:));%*(j+2)/(length(ai)+2))
wl = O2A.wl;
hold on

for k = 5:25
for j = 4:4
    %plot(wl,SRCA_SZA(:,j)*(1.02-1)/cos_zsa(j));
    plot(wl,SRCAj(:,k,j)/(ai(k)-1)/(cos_zsa(j)),'r' );
    hold on
   % plot(O2A.wl,log(O2A.piL./O2A.normpiL) - aprior*log(O2A.E./O2A.normE),'Color','r');%*(j+2)/(length(ai)+2))

    %        plot(wl,SRCAj(:,k,j)/(ai(k)-1)/cos_zsa(j));


    hold on
end
end
%set(gca, 'ylim',[0,0.1])

ylabel('residual = log(L_{norm}) - a log(E_{norm})')
xlabel('wl (nm)')
set(gca, 'FontSize',12)

%subplot(1,2,2)
%plot(O2B.wl,log(O2B.piL./O2B.normpiL) - aprior*log(O2B.E./O2B.normE),'Color',colors(b,:));%*(j+2)/(length(ai)+2))
%%set(gca, 'ylim',[0, 0.03])
%ylabel('residual = log(L_{norm}) - a log(E_{norm})')
%xlabel('wl (nm)')
%set(gca, 'FontSize',12)
%hold on
%end



%% Experiment 2. Error as function of FWHM for two values of atmospheric path length.
% this is only for O2A!
if experiment2

    for c = 1:2 % loop over two cases
        switch c
            case 1
                a = 1; % case 1, TOC with no fluorescence
            case 2
                a = 1.02; %case 2, tall tower with fluorescence
        end

        aprior = a;

        % we extrapolate the high-res spectrum to ML
        piL_MODTRAN     = piL_MODTRAN0;
        piL_MODTRAN(index2) = normpiL.* exp(log(piL_MODTRAN(index2)./normpiL)*a);

        if c>1 % add fluorescence for the second case
            piL_MODTRAN = piL_MODTRAN+ F_MODTRAN;
        end

        % convolution
        [piL,E] = deal(NaN*wl);
        [ak,Fk,FFLDk,EX]= deal(NaN*ones(length(FWHMi),1));
        for j = 1:length(FWHMi)
            FWHM = FWHMi(j);
            %FWHM = 0.3;
            sigma = FWHM/2.355; % see https://en.wikipedia.org/wiki/Full_width_at_half_maximum
            for k = 1:length(wl)
                y = normpdf(wl_MODTRAN,wl(k),sigma);
                E(k) = sum(y.*E_MODTRAN)/sum(y);
                piL(k) = sum(y.*piL_MODTRAN)/sum(y);
            end

            priorweight = 0; % no weights on the prior value (leave the retrieval free)
            %[O2A,O2B] = retrievalF(wl,E,piL,opt,aprior,0.7,1,priorweight,p,0*SRC,0*SRCB);
            [O2A,O2B] = retrievalF(wl,E,piL,opt,aprior,0.7,1,priorweight,p,0*SRCA,0*SRCB);
            Fk(j) = O2A.F;
            FFLDk(j) = O2A.iFLD;
            ak(j) = O2A.a;
            EX(j) = O2A.EXITFLAG;
        end


        figure (7)
        subplot(2,2,c),
        z=plot(FWHMi,[ak ones(length(FWHMi),1)*aprior]);
        set(z(2), 'Color','k', 'LineWidth',2)
        set(gca,'FontSize',11)
        set(gca,'ylim',[.995 1.021])
        ylabel('a')
        title(alphabet{c});

        subplot(2,2,c+2),
        z = plot(FWHMi,[Fk FFLDk ones(length(FWHMi),1)*fscope(5,760-639)*(c>1)]);
        set(z(3), 'Color','k', 'LineWidth',2)
        set(gca,'FontSize',11)
        set(gca,'ylim',[-.8 0.6])
        ylabel('F (Wm^{-2}\mum^{-1}sr^{-1})')
        %legend('Retrieved','Input')
        xlabel('FWHM (nm)')
        title(alphabet{c+2})
    end
    %subplot(313), plot(FWHMi,EX)

end

%% Experiment 3. Test this for different atmospheric profiles and elevation differences.
if experiment3
    files = {'ITALY1', 'SPAIN', 'GER1', 'GER2'};
    heights =[7,81.3; 258,324; 111,125; 069,111 ];


    for k = 1:4
        s1 = importdata(['FLEX-S3_' files{k} '.atm']);
        s2 = importdata(['FLEX-S3_' files{k} '_PC.atm']);

        for j = 1:2
            switch j
                case 1, s = s2; if k==4; s = s1; end
                case 2, s = s1; if k==4; s = s2; end
            end
            T   = s.data(:,3:20);
            atmo.M =  [T(:,1) T(:,3) T(:,4) T(:,5) T(:,12) T(:,16)];
            atmo.Ta = 25;
            [Esun,Esky] = calcIrradiance(atmo,SAIL,wl_MODTRAN);
            % figure(12)
            % plot(wl_MODTRAN,Esun), hold on
            % set(gca, 'xlim', [400 1200])

            E_MODTRAN = Esun;
            index2          = find(wl_MODTRAN>p.wl_left(1) & wl_MODTRAN<p.wl_right(1));

            % normalize the band depth by the interpolated values
            normE           = interp1([wl_MODTRAN(index2(1)),wl_MODTRAN(index2(end))],[E_MODTRAN(index2(1)),E_MODTRAN(index2(end))],wl_MODTRAN(index2));
            E_MODTRANO2A(:,j)    = E_MODTRAN(index2)./normE; %#ok<*SAGROW>
            %figure(12)
            %plot(wl_MODTRAN(index2),E_MODTRANO2A), hold on
        end
        x=log(E_MODTRANO2A(:,1));
        y=log(E_MODTRANO2A(:,2));
        slope = x \ y;
        ymod = slope*x;

        FS = interp1((640:850),fscope(5,:),wl_MODTRAN(index2));
        z = log(exp(y)+FS./normE);

        figure(12)
        subplot(2,2,k)
        plot(wl_MODTRAN(index2),[x-y,x-ymod,y-ymod]);%, z-ymod])
        xlabel('wl')
        ylabel('log(E_{norm})')
        title(alphabet{k})

        x1(k)          = (1-barometric(heights(k,2)-heights(k,1)))/barometric(heights(k,2)-heights(k,1));
        slp(k) = slope;

    end
    %
    figure(13), clf
    plot(1+x1,slp,'x','MarkerSize',8)
    hold on
    plot([1, 1.01],[1, 1.01],'k')
    xlabel('a_{barometric}')
    ylabel('a_{MODTRAN}')
end

