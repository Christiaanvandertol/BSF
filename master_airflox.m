%datapath    = 'z:\Campaign Datasets\ESA_OHP';
%clear all %#ok<CLALL>

%flox = 'ground_flox';
flox = 'airflox';
dirs        = dir( ['c:\Users\tol\projects\FRM4FLUO\' flox '\2*'] );%{'221009'};

formati     = '%2c%d%1c%d%1c%d%1c%d%1c%d%1c%d%1c';
numchar     = 4;
datetimeIDs = [7, 9, 11, 13, 15, 17]; % position of YMD-HHMMSS in FLoX


datapath    = ['c:\Users\tol\projects\FRM4FLUO\' flox];

if strcmp(flox, 'airflox')
    subfolder  = '/FLUO';
    m_name      = {'FLUO\Incoming*', 'FLUO\Reflect*'};
else
    subfolder = '';
    m_name      = {'Incoming*FLUO*', 'Reflected*FLUO*'};
end

stoptol     = 1E-9;
opt         = optimset('MaxIter',30,'TolFun',stoptol);
runSCOPE    = 0;    % this is for the comparisons to SCOPE simumlations
Lat         = 44;
Long        = 5.7;
height      = 1;
%height      = 16;%16;%60; %m
cos_vza     = 1; % set cosine of the the viewing zenith angle
x1          = (1-barometric(height))/barometric(height);
p           = set_parameters;

%load FtarbO2A.mat
%load FtarbO2B.mat
%load meanresidual2.mat
%flwf.O2A = FtarbO2A;
%flwf.O2B = FtarbO2B;

%%
load MODTRAN.mat

FLOX.FWHM = 0.31;
FLOX.wl = D(1).wl;
%[SRCA,SRCB,cO2A,cO2B] = calc_SRC(MODTRAN,FLOX,p,'test');
%[SRCA,SRCB,SRCH,cO2A,cO2B,cH2O] = calc_SRC_H2O(MODTRAN,FLOX,p,'testH2O');

load SRCA_testH2O.mat
load SRCB_testH2O.mat
load SRCH_testH2O.mat


%load SRCA_ground.mat
%load SRCB_ground.mat
%load SRCA_crane.mat
%load SRCB_crane.mat

for priorcase = 1:1
    switch priorcase
        case 1, priorweight = 0;
        case 2, priorweight = 1E12;
    end

    for k = 1:1%length(dirs)
        for m = 1:2
            files(m).fileinfo = dir([datapath '/' dirs(k).name '/' m_name{m}]);
        end
        for z = 1:length(files(1).fileinfo)
            for m = 1:2
                ifdata = 1;
                if isempty(files(1).fileinfo)
                    ifdata = 0;
                end
                if ifdata
                    for m = 1:2
                        %D(m).filename = [datapath '/' dirs{k} '/ground_flox/' fileinfo(1).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                        D(m).filename = [datapath '/' dirs(k).name subfolder '/' files(m).fileinfo(z).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                        if strcmp(flox, 'airflox')
                            [D(m).wl,D(m).data,D(m).time] = readFXBox(D(m).filename,numchar,formati,datetimeIDs);
                        else
                            [D(m).wl,D(m).data,D(m).time] = readFXBox(D(m).filename);
                        end
                    end
                    outfilename = [flox dirs(k).name];

                    %        D2(m).filename = [datapath '/' dirs{k} '/' fileinfo2(1).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                    %       [D2(m).wl,D2(m).data,D2(m).time] = readFXBox(D2(m).filename);
                end
            end
            if ifdata
                D(1).filename
                x = dirs(k).name;

                dd      = str2double(x(1:6));
                year    = 2000+floor(dd/1E4);
                month   = floor( (dd-1E4*floor(dd/1E4))/100);
                dom     = dd-1E2*floor(dd*1E-2);
                Doy     = datenum(year,month,dom)-datenum(year-1,12,31);
                cos_sza = cos(calczenithangle(Doy,24*D(m).time,0,0,Long,Lat));
                cos_sza = max(0.17,cos_sza);
                SZA     = 180/pi*acos(cos_sza);
                % the following is due to classifications of SZA of 30,45,60,75 degrees
                iSZA    = 5-max(1,min(4,(round( (min(SZA,85))/15))-1));

                E       = D(1).data;
                piL     = D(2).data;
                E(isnan(E)) = 0;
                piL(isnan(piL))=0;
                if k>10, x1 = 0; end

                if ~isempty(piL)
                    L                       = min(size(piL,2),size(E,2));
                    [Out.F]              = (nan*zeros(L,2));
                    [Out.a]              = (nan*zeros(L,2));
                    Out.time             = D(1).time;
                    %aprior                  = 1;
                    apriori                 = 1+x1*(1+cos_sza./cos_vza);
                    [Out.Efluo,Out.E0fluo]        = deal(D(1).data);
                    [Out.piLfluo,Out.piL0fluo]    = deal(D(2).data);

                    for I = 1:L
                        wl                  = D(1).wl;
                        if mean(E(:,I))>1E-3

                            % the prior information about a

                            if I>length(apriori), keyboard, end
                            aprior              = apriori(I);

                            % this is the actual retrieval
                            %keyboard
                            [O2A(I), O2B(I), H2O(I)]    = retrievalF_H2O(wl,E(:,I),piL(:,I),opt,aprior,cos_sza(I),cos_vza,priorweight,p,SRCA.*polyval(cO2A,acos(cos_sza(I))/pi*180) ,SRCB.*polyval(cO2B,acos(cos_sza(I))/pi*180),SRCH.*polyval(cH2O,acos(cos_sza(I))/pi*180));%,meanresidual); %#ok<*SAGROW>
                            %[O2A(I), O2B(I)]    = retrievalF(wl,E(:,I),piL(:,I),opt,aprior,cos_sza(I),cos_vza,priorweight,p,SRCA ,SRCB);%,
                            %Out(k).F(I,:)       = [mean(O2A(I).F) mean(O2B(I).F)];
                            Out.F(I,:)       = [(O2A(I).F) (O2B(I).F)];
                            Out.FiFLD(I,:)   = [(O2A(I).iFLD) (O2B(I).iFLD)];
                            Out.a(I,1)       = O2A(I).a;
                            Out.a(I,2)       = O2B(I).a;
                            Out.date         = dirs(k).name;
                            Out.EXITFLAG(I)  = O2A(I).EXITFLAG;
                            %aprior              = O2A(I).a;
                            %                        Out(k).wlfull       = D2(1).wl;
                            Out.wlfluo       = D(1).wl;
                            %                       Out(k).Efull        = D2(1).data;
                            %                      Out(k).piLfull      = D2(2).data
                         
                            Out.E0fluo(O2A(I).wlindex,I) = O2A(I).E0;
                            Out.piL0fluo(O2A(I).wlindex,I) = O2A(I).piL0;
                            Out.E0fluo(O2B(I).wlindex,I) = O2B(I).E0;
                            Out.piL0fluo(O2B(I).wlindex,I) = O2B(I).piL0;
                            Out.E0fluo(H2O(I).wlindex,I) = H2O(I).E0;
                            Out.piL0fluo(H2O(I).wlindex,I) = H2O(I).piL0;
                            
                            Out.O2A(I)       = O2A(I);
                            Out.O2B(I)       = O2B(I);
                            Out.H2O(I)       = H2O(I);


                            Out.cos_sza      = cos_sza;
                        end
                    end
                end

                % save the output

                save(['../output/' num2str(priorcase) outfilename '.mat'], 'Out')
                csvoutput = [Out.time*24, Out.F];
                csvwrite(['../output/' num2str(priorcase) outfilename '.csv'],csvoutput)
                %clear('Out')
            end
        end


    end

end

%%




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
%             save('c:\Users\tol\Documents\models\retrieval_develop\data\measured\ATMOFLEX\r.txt','r','-ascii')
%             save('c:\Users\tol\Documents\models\retrieval_develop\data\measured\ATMOFLEX\wl.txt','wl','-ascii')
%
%             %% run retrieval
%             % retrieve the vegetation properties from the measured reflectance
%             run('c:\Users\tol\Documents\models\retrieval_develop\canopy_soil retrieval\master')
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
%         rad.Esun_ = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\Esun.csv'],',',2,0);
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








