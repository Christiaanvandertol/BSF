%%
wl_left     = [759, 686.5]; % just left, used in new method
wl_right    = [768, 695.0]; % just right, used in new method
stoptol     = 1E-4;
opt         = optimset('MaxIter',30,'TolFun',stoptol);

clear input
figure(10), clf
%J = 494;
J = 400;
Day = 4;
for O2band = 1:1
    switch O2band
        case 1
            data = Out(Day).O2A(J);
        case 2
            data = Out(Day).O2B(J);
            
    end
    
    wl = data.wl;
    E = data.E;
    piL = data.piL;
    piLr = data.piLr;
    normE = data.normE;
    normpiL = data.normpiL;
    ar = Out(Day).a(J);
    index2          = find(wl>wl_left(O2band) & wl<wl_right(O2band));

    switch O2band
        case 1
            input.flwf      = 0.66+0.34*(length(index2):-1:1)'/length(index2);
            input.atcor     = 1;
            input.logxlim   = -.7;
        case 2
            input.flwf      = 1;
            input.atcor     = 1;
            input.logxlim   = -.3;
           % input.a         = ar2;
    end
    
   % normE       	= interp1([wl(index2(1)),wl(index2(end))],[E(index2(1)),E(index2(end))],wl(index2));
  %  normpiL         = interp1([wl(index2(1)),wl(index2(end))],[piL(index2(1)),piL(index2(end))],wl(index2));
    
    % normalize the band depth by the interpolated values
    input.logx      = log(E(index2)./normE);%+ logtaua;
    input.y         = piL(index2)./normpiL;%*taua;
    input.normpiL 	= normpiL;
    %input.aprior    = aprior;
    input.aprior    = Out(Day).a(J-1);
    input.priorweight = 0;
    input.cos_sza   = Out(Day).cos_sza(J);
    input.cos_vza   = cos_vza;
    input.SRC       = SRC;
    
    Fj = (0:.01:1)' *1E-3;
    am = (1:.001:1.04)';
    RMSE = Fj*NaN;
    [Nb,a] = deal(length(Fj));
    Na = length(am);
    for k = 1:Nb
        [Er,a(k)] = cost4F(Fj(k),input);
        RMSE(k) = sqrt((Er'*Er)/length(Er));
    end
    
    input.atcor = 0;
    RMSEi = zeros(Nb,Na);
    for k = 1:Nb
        for m = 1:Na
            input.a = am(m);
            [Er] = cost4F(Fj(k),input);
            RMSEi(k,m) = sqrt((Er'*Er)/length(Er));
        end
    end
    
        
    
    f           = @(Fi)cost4F(Fi,input);
    [Fr,RESNORM,RESIDUAL,EXITFLAG]           = lsqnonlin(f,0,-.1,.1,opt); %#ok<ASGLU> % F normalized by the radiance outside the band
    [Erm,ar2,piLr] = cost4F(Fr,input);
    
    
    
    %subplot(2,3,O2band*3-2)
    subplot(312)
    plot(Fj*1E3,RMSE,'k')
    xlabel('F_{760} (Wm^{-2}\mum{-1}sr^{-1})')
    ylabel('RMSE')
    
    
    %subplot(2,3,O2band*3-1)
    subplot(313)
    plot(Fj*1E3,a,'k')
    hold on
    [~,I] = min(abs(Fj-(data.F-0.2E-3)));
    [~,J] = min(abs(Fj-(data.F+0.2E-3)));
    [~,K] = min(abs(Fj-(data.F)));
    plot([data.F*1E3-0.2,data.F*1E3-0.2],[1,a(I)],'k--')
    plot([0,data.F*1E3-0.2],[a(I),a(I)],'k--')
    plot([data.F*1E3+0.2,data.F*1E3+0.2],[1,a(J)],'k--')
    plot([0,data.F*1E3+0.2],[a(J),a(J)],'k--')
    plot([data.F*1E3,data.F*1E3],[1,a(K)],'r--')
    plot([0,data.F*1E3],[a(K),a(K)],'r--')
    
    xlabel('F_{760} (Wm^{-2}\mum{-1}sr^{-1})')
    ylabel('a')
    
%     %subplot(2,3,O2band*3)
%     subplot(224)
%     plot(log(E./normE), log(piL./normpiL), 'kx','MarkerSize',3)
%     hold on
%     plot([-2.5 0],[-2.5 0], 'k--')
%     %   if O2band ==1
%     plot(log(E./normE), log(piLr), 'r*', 'MarkerSize',3)
%     plot([-2.5 0],ar2*[-2.5 0],'r')
%     if O2band ==2, set(gca,'xlim', [-1 0],'ylim', [-1 0]),
%     else
%         set(gca,'xlim', [-3 0],'ylim', [-3 0]),
%     end
%     xlabel('log E(\lambda)/log(E_0)')
%     ylabel('log \piL(\lambda)/log(rE_0)')
    
    subplot(311)
    f = contourf(Fj*1E3*ones(1,Na),ones(Nb,1)*am',RMSEi,'LineColor', 'none','LevelList',-6:0.001: 0.035);
    hold on
    plot(Fj*1E3,a','k')
    hold on
    %plot([Fr*1E3+.2 Fr*1E3+.2],[1 ar],'k')
    xlabel('F_{760} (Wm^{-2}\mum{-1}sr^{-1})')
    ylabel('a')
    colorbar
    
end

%%

