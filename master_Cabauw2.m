%datapath    = 'z:\Campaign Datasets\ESA_OHP';
%clear all %#ok<CLALL>

%dirs        = dir('c:\Users\tol\Documents\projects\DEFLOX\ESA_DEFLOX_CCN4\DATA\cabauw\flox016m\2*' );%{'221009'};
dirs        = dir('c:\Users\tol\projects\FRM4FLUO\Cabauw\flox100m\2*' );%{'221009'};
%dirs        = dir('y:\projects\ESA_DEFLOX_CCN4\DATA\cabauw\flox016m\2*' );%{'221009'};

for d =1:length(dirs)
    for flx = 2:2
        switch flx
            case 1, flox        = 'flox200m'; height = 200; clear('Out')
            case 2, flox        = 'flox100m'; height = 100; clear('Out')
            case 3, flox        = 'flox016m'; height = 16; clear('Out')
        end

        %flox        = 'ground_Flox
        %datapath    = 'y:\projects\ESA_DEFLOX_CCN4\Data\Cabauw\';
        datapath    = 'c:\Users\tol\projects\FRM4FLUO\Cabauw\';
        m_name      = {'Incoming*fluo*', 'Reflect*fluo*','Incoming*full*', 'Reflect*full*'};
        stoptol     = 1E-4;
        stopx       = 1E-3;
        opt         = optimset('MaxIter',30,'TolFun',stoptol,'TolX',stopx);
        runSCOPE    = 0;    % this is for the comparisons to SCOPE simumlations
        Lat         = 52.0;
        Long        = 4.9;
        %height      = 16;%16;%60; %m
        cos_vza     = cos(13/180*pi);%1; % set cosine of the the viewing zenith angle
        x1          = (1-barometric(height))/barometric(height);
        p           = set_parameters;

        %% tarb fluorescence spectrum

        %    load FtarbO2A.mat
        %    load FtarbO2B.mat
        %   flwf.O2A = FtarbO2A;
        %  flwf.O2B = FtarbO2B;

        %%
        load MODTRAN.mat

        FLOX.FWHM = 0.31;
        FLOX.wl = D(1).wl;
        [SRCA,SRCB,cO2A,cO2B] = calc_SRC(MODTRAN,FLOX,p,'cabauw100');
    %    [SRCA,SRCB,SRCH,cO2A,cO2B,cH2O] = calc_SRC_H2O(MODTRAN,FLOX,p,'cabauw100');

        %load('y:\projects\ESA_DEFLOX_CCN4\BSF\SRCA_ground.mat')
        load('SRCA_cabauw100.mat')
        load('SRCB_cabauw100.mat')
        %cO2A = c;
        %load SRCA_BRDF.mat
        %load('y:\projects\ESA_DEFLOX_CCN4\BSF\SRCB_ground.mat')
        %cO2B = c;

        %  SCRA = 0*SRCA;

        for priorcase = 1:1
            switch priorcase
                case 1, priorweight = 0;
                case 2, priorweight = 1E12;
            end

            %        for k = 1:length(dirs)
            for m = 1:4
                files(m).fileinfo = dir([datapath flox '/' dirs(d).name '/' m_name{m}]);
            end
            for z = 1:length(files(1).fileinfo)
                for m = 1:2
                    ifdata = 1;
                    if isempty(files(1).fileinfo)
                        ifdata = 0;
                    end
                    if ifdata
                        for m = 1:4
                            %D(m).filename = [datapath '/' dirs{k} '/ground_flox/' fileinfo(1).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                            D(m).filename = [datapath flox '/' dirs(d).name '/' files(m).fileinfo(z).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                            [D(m).wl,D(m).data,D(m).time] = readFXBox(D(m).filename,4+2*(m>2));
                            % keyboard
                        end
                        outfilename = ['cabauw' flox dirs(d).name];

                        %        D2(m).filename = [datapath '/' dirs{k} '/' fileinfo2(1).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                        %       [D2(m).wl,D2(m).data,D2(m).time] = readFXBox(D2(m).filename);
                    end
                end
                if ifdata
                    D(1).filename
                    x = dirs(d).name;
                    Ei      = D(1).data;
                    piLi    = D(2).data;

                    Ei      = Ei(:,1:length(D(1).time));
                    piLi    = piLi(:,1:length(D(2).time));

                    % filter on NDVI, to remove some erroneous FloX data
                    %            C       = find(Ei(600,:)./Ei(150,:)<1);

                    %           E       = Ei(:,C);
                    %           piL     = piLi(:,C);
                    %                    if k>10, x1 = 0; end

                    E = Ei; % comment out for NDVI filtering
                    piL = piLi;

                    dd      = str2double(x(1:6));
                    year    = 2000+floor(dd/1E4);
                    month   = floor( (dd-1E4*floor(dd/1E4))/100);
                    dom     = dd-1E2*floor(dd*1E-2);
                    Doy     = datenum(year,month,dom)-datenum(year-1,12,31);
                    %                    cos_sza = cos(calczenithangle(Doy,24*D(m).time(C),0,0,Long,Lat));
                    cos_sza = cos(calczenithangle(Doy,24*D(1).time,0,0,Long,Lat));
                    cos_sza = max(0.17,cos_sza);
                    SZA     = 180/pi*acos(cos_sza);
                    % the following is due to classifications of SZA of 30,45,60,75 degrees
                    iSZA    = 5-max(1,min(4,(round( (min(SZA,85))/15))-1));




                    if ~isempty(piL)
                        L                       = min(size(piL,2),size(E,2));
                        [Out.F]              = (nan*zeros(L,2));
                        [Out.a]              = (nan*zeros(L,2));
                        %Out.time             = D(1).time(C);
                        Out.time             = D(1).time;
                        %apriori                  = 1+x1*(1+cos_sza./cos_vza)*0;
                        apriori                 = 1+x1*(1+cos_sza./cos_vza);
                     z%   aprior = 0*apriori+1.03;
                                            [Out.Efluo,Out.E0fluo]        = deal(D(1).data);
                    [Out.piLfluo,Out.piL0fluo]    = deal(D(2).data);

                        for I = 1:L%129:136%L
                            wl                  = D(1).wl;
                            if mean(E(:,I))>1E-3

                                % the prior information about a

                                %        if I>length(apriori), keyboard, end
                                aprior              = apriori(I);

                                % this is the actual retrieval
                                %keyboard
                                J = max(1,I-5):min(L,I+5);
                                J = I;%
                                %[O2A(I), O2B(I)]    = retrievalF(wl,mean(E(:,J),2),mean(piL(:,J),2),opt,aprior,cos_sza(I),cos_vza,priorweight,p,SRCA/cos_sza(I),SRCB/cos_sza(I));%,meanresidual); %#ok<*SAGROW>
                                %keyboard

                                [O2A(I), O2B(I)]    = retrievalF(wl,mean(E(:,J),2),mean(piL(:,J),2),opt,aprior,cos_sza(I),cos_vza,priorweight,p,SRCA.*polyval(cO2A,acos(cos_sza(I))/pi*180) ,SRCB.*polyval(cO2B,acos(cos_sza(I))/pi*180) );%,meanresidual); %#ok<*SAGROW>
                                %                          keyboard
                                %[O2A(I), O2B(I)]    = retrievalF(wl,E(:,I),piL(:,I),opt,aprior,cos_sza(I),cos_vza,priorweight,p,SRCA/cos_sza(I),SRCB/cos_sza(I));%,meanresidual); %#ok<*SAGROW>

                                %Out(k).F(I,:)       = [mean(O2A(I).F) mean(O2B(I).F)];

                                Out.F(I,:)       = [(O2A(I).F) (O2B(I).F)];
                                Out.FiFLD(I,:)   = [(O2A(I).iFLD) (O2B(I).iFLD)];
                                Out.a(I,1)       = O2A(I).a;
                                Out.a(I,2)       = O2B(I).a;
                                Out.date         = dirs(d).name;
                                Out.EXITFLAG(I)  = O2A(I).EXITFLAG;
                                Out.RESIDUAL(I)  = mean(abs(O2A(I).RESIDUAL));
                                %aprior              = O2A(I).a;
                                %                        Out(k).wlfull       = D2(1).wl;
                                Out.wlfluo       = D(1).wl;
                                %                       Out(k).Efull        = D2(1).data;
                                %                      Out(k).piLfull      = D2(2).data;
                                Out.Efluo        = E;%D(1).data;
                                Out.piLfluo      = piL;%D(2).data;
                                Out.O2A(I)       = O2A(I);
                                Out.O2B(I)       = O2B(I);
                                Out.cos_sza      = cos_sza;
                                Out.wlfull       = D(3).wl;
                                Out.Efull        = D(3).data;
                                Out.piLfull      = D(4).data;
                                Out.timefull     = D(3).time;

                                                            Out.E0fluo(O2A(I).wlindex,I) = O2A(I).E0;
                            Out.piL0fluo(O2A(I).wlindex,I) = O2A(I).piL0;
                            Out.E0fluo(O2B(I).wlindex,I) = O2B(I).E0;
                            Out.piL0fluo(O2B(I).wlindex,I) = O2B(I).piL0;
                                %                      keyboard
                            end
                        end
                    end
                    % save the output

                    % save(['../output/prior_noSRA' num2str(priorcase) outfilename '.mat'], 'Out')
                    %  csvoutput = [Out.time*24, Out.F];
                    %  csvwrite(['../output/prior_noSRA' num2str(priorcase) outfilename '.csv'],csvoutput)
                    %clear('Out')
                end
            end
        end
    end
    Out_all(d) = Out;
end


%         %% Comparison against SCOPE output
%
%         if runSCOPE
%             r = zeros(length(Out(k).piLfull),length(dirs));
%             for k = 1:length(dirs)
%                 r(:,k) = mean(Out.piLfull,2)./mean(Out.Efull,2);
%             end
%             if k == 10 % this was a different sprectrometer (different Floxbox!)
%                 r(:,10) = interp1(Out.wlfull, r(:,10),Out.wlfull);
%             end
%
%             wl = Out.wlfull;
%             % save('c:\Users\tol\Documents\models\retrieval_develop\data\measured\ATMOFLEX\r.txt','r','-ascii')
%             % save('c:\Users\tol\Documents\models\retrieval_develop\data\measured\ATMOFLEX\wl.txt','wl','-ascii')

%             save('c:\Users\tol\retrieval_RTMO-master\data\measured\ATMOFLEX\r.txt','r','-ascii')
%             save('c:\Users\tol\retrieval_RTMO-master\models\retrieval_develop\data\measured\ATMOFLEX\wl.txt','wl','-ascii')
%
%             %% run retrieval
%             % retrieve the vegetation properties from the measured reflectance
%             % run('c:\Users\tol\Documents\models\retrieval_develop\canopy_soil retrieval\master')
%             %% run SCOPE
%             % place the retrieved values in the SCOPE input
%             run('..\code\SCOPE2\SCOPE')
%         end
%
%         %     % after running SCOPE, specify here the output folder
%         SCOPEoutputfolder = 'OHP_2021-06-25-1444';
%         %    SCOPEoutputfolder = 'OHP_2022-07-10-2357';
%
%         fscope = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\fluorescence.csv'],',',2,0);
%         spectral.wlS = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\wlS.txt'],',',0,0);
% %         Esun_ = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\Esun.csv'],',',2,0);
%         rad.Esky_ = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\Esky.csv'],',',2,0);
%         rad.Esun_ = rad.Esun_';
%         rad.Esky_ = rad.Esky_';
%
%         %%
%         for k = 1:length(dirs)
%
%             % scaling of the output of SCOPE to the measured irradiance
%             % (described in the paper)
%             %        iwlQ = find(Out(k).wlfull>400 & Out(k).wlfull<900);
%             iwlQs = find(spectral.wlS>400 & spectral.wlS<900);
%             %        wlQ  = Out(k).wlfull(iwlQ);
%             %       EQ = Sint(Out(k).Efull(iwlQ,:)',wlQ);
%             %        EQS = 1E-3*Sint(rad.Esun_(iwlQs,1)+rad.Esky_(iwlQs,1),spectral.wlS(iwlQs));
%             %       M = pi*EQ./EQS;
%             M = 1;
%             if k<11
%                 Out.Fscope = [fscope(k,759-639)*M,mean(fscope(k,685-639:700-639))*M];
%             elseif k == 11 % this is the same day as day 5 (see 'dirs' in one of the first lines)
%                 Out.Fscope = [fscope(5,759-639)*M,mean(fscope(5,685-639:700-639))*M];
%             else
%                 Out.Fscope = [fscope(10,759-639)*M,mean(fscope(10,685-639:700-639))*M];
%             end








