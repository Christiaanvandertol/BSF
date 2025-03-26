function [SRCA,SRCB,cO2A,cO2B] = calc_SRC(MODTRAN,FLOX,p,outfilename)

%% parameters for the retrieval algorithm
stoptol     = 1E-6;
opt         = optimset('MaxIter',30,'TolFun',stoptol);

%%
wl          = FLOX.wl;
FWHM        = FLOX.FWHM;
SZA         = MODTRAN.SZA;
%cos_sza     = cos(SZA/180*pi);

%% loading SCOPE reflectances and fluorescence
%SCOPEoutputfolder = 'OHP_2021-06-25-1444';
%fscope      = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\fluorescence.csv'],',',2,0);
%wlS         = load(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\wlS.txt']);


%% For an input value of FWMH, error as function of atmospheric path length
% In this piece of code, we also estimate a correction function for the
% effect of the SRF on the spectral shape of the O2A and O2B bands. This
% is estimated as follows:
% We estimate <L> as:
% Second: we estimate L_{approx} as: exp(log(<E>r)*a
% Third: SRC = (<L>-L_{approx}) / (1-a)

%ai          = (1.000:.004:1.1); % this was before 11 March 2025
ai          = (1.000:.04:1.5); % changed this on 11 March 2025
priorweight = 0;
sigma       = FWHM/2.355; % see https://en.wikipedia.org/wiki/Full_width_at_half_maximum

%% initialization
for O2band = 1:2
    O2(O2band).index  = find(wl>p.wl_left(O2band) & wl<p.wl_right(O2band)); %#ok<*AGROW>
    O2(O2band).SRC_SZA     = zeros(length(O2(O2band).index),4);
    O2(O2band).SRCj        = zeros(length(O2(O2band).index),length(ai),4);
    O2(O2band).SRC         = zeros(length(O2(O2band).index),1);
end
[piL,E]     = deal(NaN*wl);
[ak,Fk,FFLDk,EX]= deal(NaN*ones(length(ai),1));

%%
for iSZA = 1:length(SZA)
    wl_MODTRAN      = MODTRAN.wl(:,iSZA);
    piL_MODTRAN0    = MODTRAN.piL(:,iSZA);
    E_MODTRAN       = MODTRAN.E(:,iSZA);    
    normpiL2        = MODTRAN.normpiO2A(:,iSZA);
    normpiL3        = MODTRAN.normpiO2B(:,iSZA);
    index2          = MODTRAN.iO2A(:,iSZA);
    index3          = MODTRAN.iO2B(:,iSZA);

    for j = 1:length(ai)
        a                   = ai(j);
        aprior              = a;
        piL_MODTRAN         = piL_MODTRAN0;
        piL_MODTRAN(index2) = normpiL2.* exp(log(piL_MODTRAN(index2)./normpiL2)*a);
        piL_MODTRAN(index3) = normpiL3.* exp(log(piL_MODTRAN(index3)./normpiL3)*a);

        % convolution
        for k = 1:length(wl)
            y               = normpdf(wl_MODTRAN,wl(k),sigma);
            E(k)            = sum(y.*E_MODTRAN)/sum(y);
            piL(k)          = sum(y.*piL_MODTRAN)/sum(y);
        end

        %[O2A,O2B] = retrievalF(wl,E,piL,opt,aprior,0.7,1,priorweight,p,O2(1).SRC*0,O2(2).SRC*0);
        [O2A,O2B] = retrievalF(wl,E,piL,opt,aprior,cos(SZA(iSZA)/180*pi),1,priorweight,p,O2(1).SRC*0,O2(2).SRC*0);

        Fk(j)       = O2A.F;
        FFLDk(j)    = O2A.iFLD;
        ak(j)       = O2A.a;
        EX(j)       = O2A.EXITFLAG;

        O2(1).SRCj(:,j,iSZA) = (log(O2A.piL./O2A.normpiL) - aprior* log(O2A.E./O2A.normE));
        O2(2).SRCj(:,j,iSZA) = (log(O2B.piL./O2B.normpiL) - aprior* log(O2B.E./O2B.normE));
    end
    for O2band = 1:2
        for k = 1:size(O2(O2band).SRC,1)
            O2(O2band).SRC_SZA(k,iSZA) = (ai-1)' \ O2(O2band).SRCj(k,:,iSZA)';
            %alternative:                c = polyfit(ai-1,SRCBj(k,:,iSZA),1); alternative to linear
            %                regression. With an offset, but this is not needed...
        end
    end
end


%SRCA = SRCA_SZA(:,4)./cos(30/180*pi).^4;
for O2band = 1:2
    O2(O2band).SRC  = O2(O2band).SRC_SZA(:,4)./mean(O2(O2band).SRC_SZA(:,4));
    O2(O2band).c    = polyfit((SZA), mean(O2(O2band).SRC_SZA),2);
    %O2(O2band).c    = polyfit((SZA), mean(O2(O2band).SRC_SZA),1);
    %keyboard
end

SRCA = O2(1).SRC;
SRCB = O2(2).SRC;
cO2A = O2(1).c;
cO2B = O2(2).c;

save(['SRCA_' outfilename '.mat'],'SRCA', 'cO2A')
save(['SRCB_' outfilename '.mat'],'SRCB', 'cO2B')
