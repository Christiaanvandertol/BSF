%datapath    = 'z:\Campaign Datasets\ESA_OHP';

datapath    = '..\data\OHP\';
dirs      = {'180319', '180420','180425','180525','180927','190321','190713','190902','200624','180927b','180927b','180927'};
m_name      = {'Incoming', 'Reflect*'};
stoptol     = 1E-4;
opt         = optimset('MaxIter',30,'TolFun',stoptol);
runSCOPE    = 0;    % this is for the comparisons to SCOPE simumlations
Lat         = 44;
Long        = 5.7;
height      = 100; %m
cos_vza     = 1; % set cosine of the the viewing zenith angle
x1          = (1-barometric(height))/barometric(height);
p           = set_parameters;
load SRCA.mat
load SRCB.mat

for iprior = 1:1
    switch iprior
        case 1
            priorweight = 0;
            outfilename = 'Noprior.mat';
        case 2
            priorweight = 1E0;
            outfilename = 'Prior.mat';
        case 3
            priorweight = 1E12;
            outfilename = 'Forced.mat';
    end
    
    for k = 1:length(dirs)
        for m = 1:2
            fileinfo = dir([datapath '/' dirs{k} '/' m_name{m} '_radiance_FLUO_*.csv']);
            fileinfo2 = dir([datapath '/' dirs{k} '/' m_name{m} '_radiance_FULL_*.csv']);
            
            ifdata = 1;
            if isempty(fileinfo)
                ifdata = 0;
            end
            if ifdata
                D(m).filename = [datapath '/' dirs{k} '/' fileinfo(1).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                [D(m).wl,D(m).data,D(m).time] = readFXBox(D(m).filename);
                D2(m).filename = [datapath '/' dirs{k} '/' fileinfo2(1).name];%#ok<*SAGROW> %[datapath '/' dirs(k).name '/' m_name{m} '_radiance_FLUO_*.csv']; %#ok<*SAGROW>
                [D2(m).wl,D2(m).data,D2(m).time] = readFXBox(D2(m).filename);
            end
        end
        if ifdata
            D(1).filename
            x = dirs{k};
            
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
            if k>10, x1 = 0; end
            
            if ~isempty(piL)
                L                       = min(size(piL,2),size(E,2))-1;
                [Out(k).F]              = (nan*zeros(L,2));
                [Out(k).a]              = (nan*zeros(L,2));
                Out(k).time             = D(1).time;
                %aprior                  = 1;
                apriori                 = 1+x1*(1+cos_sza./cos_vza);
                for I = 1:L
                    wl                  = D(1).wl;
                    if mean(E(:,I))>1E-3
                        
                        % the prior information about a
                        
                        aprior              = apriori(I);

                        % this is the actual retrieval
                        [O2A(I), O2B(I)]    = retrievalF(wl,E(:,I),piL(:,I),opt,aprior,cos_sza(I),cos_vza,priorweight,p,SRCA/cos_sza(I),SRCB/cos_sza(I)); %#ok<*SAGROW>
                        %Out(k).F(I,:)       = [mean(O2A(I).F) mean(O2B(I).F)];
                        Out(k).F(I,:)       = [(O2A(I).F) (O2B(I).F)];
                        Out(k).FiFLD(I,:)   = [(O2A(I).iFLD) (O2B(I).iFLD)];
                        Out(k).a(I,1)       = O2A(I).a;
                        Out(k).a(I,2)       = O2B(I).a;
                        Out(k).date         = dirs{k};
                        Out(k).EXITFLAG(I)  = O2A(I).EXITFLAG;
                        %aprior              = O2A(I).a;
                        Out(k).wlfull       = D2(1).wl;
                        Out(k).wlfluo       = D(1).wl;
                        Out(k).Efull        = D2(1).data;
                        Out(k).piLfull      = D2(2).data;
                        Out(k).Efluo        = D(1).data;
                        Out(k).piLfluo      = D(2).data;
                        Out(k).O2A(I)       = O2A(I);
                        Out(k).O2B(I)       = O2B(I);
                        Out(k).cos_sza      = cos_sza;
                    end
                end
            end
        end
    end
    
    %% Comparison against SCOPE output
    
    if runSCOPE
        r = zeros(length(Out(k).piLfull),length(dirs));
        for k = 1:length(dirs)
            r(:,k) = mean(Out(k).piLfull,2)./mean(Out(k).Efull,2);
        end
        if k == 10 % this was a different sprectrometer (different Floxbox!)
            r(:,10) = interp1(Out(10).wlfull, r(:,10),Out(1).wlfull);
        end
        
        wl = Out.wlfull;
        save('c:\Users\tol\Documents\models\retrieval_develop\data\measured\ATMOFLEX\r.txt','r','-ascii')
        save('c:\Users\tol\Documents\models\retrieval_develop\data\measured\ATMOFLEX\wl.txt','wl','-ascii')
        
        %% run retrieval
        % retrieve the vegetation properties from the measured reflectance
        run('c:\Users\tol\Documents\models\retrieval_develop\canopy_soil retrieval\master')
        %% run SCOPE
        % place the retrieved values in the SCOPE input
        run('..\code\SCOPE2\SCOPE')
    end
    
    % after running SCOPE, specify here the output folder
    SCOPEoutputfolder = 'OHP_2021-06-25-1444';
%    SCOPEoutputfolder = 'OHP_2022-07-10-2357';
    
    fscope = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\fluorescence.csv'],',',2,0);
    spectral.wlS = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\wlS.txt'],',',0,0);
    rad.Esun_ = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\Esun.csv'],',',2,0);
    rad.Esky_ = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\Esky.csv'],',',2,0);
    rad.Esun_ = rad.Esun_';
    rad.Esky_ = rad.Esky_';
    
    %%
    for k = 1:length(dirs)
        
        % scaling of the output of SCOPE to the measured irradiance
        % (described in the paper)
        iwlQ = find(Out(k).wlfull>400 & Out(k).wlfull<900);
        iwlQs = find(spectral.wlS>400 & spectral.wlS<900);
        wlQ  = Out(k).wlfull(iwlQ);
        EQ = Sint(Out(k).Efull(iwlQ,:)',wlQ);
        EQS = 1E-3*Sint(rad.Esun_(iwlQs,1)+rad.Esky_(iwlQs,1),spectral.wlS(iwlQs));
        M = pi*EQ./EQS;
        if k<11
            Out(k).Fscope = [fscope(k,759-639)*M,mean(fscope(k,685-639:700-639))*M];
        elseif k == 11 % this is the same day as day 5 (see 'dirs' in one of the first lines)
            Out(k).Fscope = [fscope(5,759-639)*M,mean(fscope(5,685-639:700-639))*M];
        else
            Out(k).Fscope = [fscope(10,759-639)*M,mean(fscope(10,685-639:700-639))*M];
        end
    end
    
    % save the output
    save(['../output/' outfilename], 'Out');
    
end




