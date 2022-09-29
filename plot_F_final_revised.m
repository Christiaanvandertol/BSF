dates = {...
    '19 March 2018' ...
    '20 April 2018' ...
    '25 April 2018' ...
    '25 May 2018' ...
    '27 Sept 2018'...
    '21 March 2019' ...
    '13 July 2019' ...
    '2 Sept 2019' ...
    '24 June 2020' ...
    '27 Sept 2018' ...
    '27 Sept 2018'};

O2bandname = {'O2A', 'O2B'};
MS = 2;
Colors = {'b', 'r', 'm'};
ColorWhite = [.8 .8 1; 1 .8 .8; 1 .6 1];
Markers = {'x','o','v'};
alphabet = char({'abcdefghikjlmnopqrstuvwxyz'});

cos_vza     = 1;
height      = 100; %m
x1 = (1-barometric(height))/barometric(height);
constants=define_constants;
load('../output/Noprior.mat'); Noprior = Out;
load('../output/Prior.mat'); Weakprior = Out;
load('../output/Forced.mat'); Forced = Out;

m=1;
for n = 1:2
    Out = Noprior;

    O2band = n; % added 23/9/2022 CvdT
    %for m = 1:3
%    switch m
%        case 1, Out = Forced;
%        case 3, Out = Weakprior;
%        case 2, Out = Noprior;
%    end
    Color = Colors{m};
    Marker = Markers{m};
    
    for d = 1:9%length(dates)
        Datestr = datestr(dates{d},'yymmdd');
        
        F = Out(d).F;
        a = Out(d).a;
        F(imag(F)>0) = NaN;
        a(imag(a)>0)= NaN;
        
        F1 = figure(1);
        set(F1, 'Position', [488.2000 175.4000 739.2000 586.4000]);
        subplot(3,3,d)
        hold on
        if n == 1 
            plot(Out(d).time(10:end-11)*24,a(10:end-10,1),'o','MarkerSize', MS,'MarkerEdgeColor',ColorWhite(m,:)), hold on
        else
            plot(Out(d).time(10:end-11)*24,movmean(a(10:end-10,1),20),'b-','Color',Color,'LineWidth',1);
            
        end
        plot(Out(d).time(3:end-1)*24, 1 + x1*(1+Out(d).cos_sza(3:end-1)./cos_vza),'k')
        set(gca, 'ylim', [ 1 1.05], 'xlim', [6 18], 'xtick', 0:6:24)
        if d>6, xlabel('hour'), end
        if mod(d,3)==1, ylabel('a'), end
        title([dates{d} ' ' O2bandname(O2band)])
        
        for O2band = 1:2
            Fn= figure(O2band+1) ;
            set(Fn, 'Position', [488.2000 175.4000 739.2000 586.4000]);
            subplot(3,3,d)
            hold on
            if n == 1
                plot(Out(d).time(10:end-11)*24,F(10:end-10,O2band)*1E3,'bo', 'MarkerSize', MS,'MarkerEdgeColor',ColorWhite(m,:)); hold on
            else
                plot(Out(d).time(10:end-11)*24,movmean(F(10:end-10,O2band)*1E3,20),'b-','Color',Color,'LineWidth',1);
                %if m==2, keyboard, end
            end
            plot(Out(d).time(2:end-1)*24,Out(d).Fscope(2:end-1,O2band),'LineWidth',1,'Color', 'k');
            
            set(gca, 'ylim', [ -.4 1.4], 'xlim', [6 18], 'xtick', 0:6:24)
            title(dates{d})
            if d>6, xlabel('hour'), end
            if mod(d,3)==1, ylabel('F (Wm^{-2}\mum^{-1}sr^{-1})'), end
            %     if d == 1, legend(['retrieved ' char(O2bandname(O2band))],['SCOPE '  char(O2bandname(O2band))] ), end
        end
  
        F7 = figure(7) ;
        set(F7, 'Position', [488.2000 175.4000 739.2000 586.4000]);
        subplot(3,3,d)
        mn = mean(Out(d).piLfull(:,200:end-250),2)./mean(Out(d).Efull(:,200:end-250),2);
        plot(Out(d).wlfull, mn,'k')
        hold on
        st = std( Out(d).piLfull(:,200:end-250)./Out(d).Efull(:,200:end-250),[],2);
        
        plot(Out(d).wlfull, mn+st,'color',[.4 .4 .4])
        plot(Out(d).wlfull, mn-st,'color',[.4 .4 .4])
        set(gca, 'ylim',[0,1])
        %plot(Out(k).wlfull, (Out(k).piLfull)./(Out(k).Efull))
        if d>6
            xlabel('wl (nm)')
        end
        if mod(d,3)==1
            ylabel('r')
        end
        
        set(gca,'xlim', [400,900])
        title(dates{d})
        set(gca, 'ylim', [0,.5])
        
        F8 = figure(8) ;
        set(F8, 'Position', [488.2000 175.4000 739.2000 586.4000]);
        subplot(3,3,d)
        I = find(Out(d).wlfull>400 & Out(d).wlfull<700);
        PAR = zeros(size(Out(d).Efull,2),1);
        for k = 1:size(Out(d).Efull,2)
            PAR(k) = 1E6*Sint(Out(d).wlfull(I),e2phot(1E-9*Out(d).wlfull(I),pi*Out(d).Efull(I,k),constants));
        end
        plot(Out(d).time(1:end-1)*24, PAR(1:end-1),'k')
        %plot(Out(k).wlfull, (Out(k).piLfull)./(Out(k).Efull))
        % if d>6
        %    xlabel('wl (nm)')
        %end
        %    if mod(d,3)==1
        %       ylabel('r')
        %  end
        
        % set(gca,'xlim', [400,900])
        set(gca, 'ylim', [ 0 2000], 'xlim', [6 18], 'xtick', 0:6:24)
        title(dates{d})
        if d>6, xlabel('hour'), end
        if (mod(d,3)==1), ylabel('iPAR (\mumol m^{-2}s^{-1}'), end
        %set(gca, 'ylim', [0,700])
    end
end
%end

%%
%Out = Noprior;


d = [11,12,11,5];
titles = {'low / no atmos char', 'high / no atmos char','low / atmos char','high / atmos char'};

h = [0,x1, 0, x1];

for k = 1:4
    if k<3
        Out = Forced(d(k));
    else
        Out = Noprior(d(k));
    end
    
    F = Out.F;
    a = Out.a;
    F(imag(F)>0)    = NaN;
    a(imag(a)>0)    = NaN;
    
    for O2band = 1:2
        F4 = figure(14+O2band);
        switch O2band
            case 1, set(F4, 'Position', [488.2000 189.8000 508.8000 572]);
            case 2, set(F4, 'Position', [488.2000 384.2000 507.2000 377.6000]);
        end
        
        if k>2 && O2band==1
            subplot(3,2,k+2)
            hold on
            plot(Out.time(3:end-1)*24,a(3:end,1),'bx','Marker', Marker,'MarkerSize', MS,'MarkerEdgeColor',[.5 .5 1]), hold on% ColorWhite(:,1)), hold on
            plot(Out.time(10:end-11)*24,movmean(a(10:end-10,1),20),'b','LineWidth',1), hold on
            plot(Out.time(3:end-1)*24, 1 + h(k)*(1+Out.cos_sza(3:end-1)./cos_vza),'k')
            set(gca, 'ylim', [ 0.96 1.05], 'xlim', [6 18], 'xtick', 0:6:24)
            xlabel('hour')
            ylabel('a')
            

            % title(titles{k})
        end
        
        %    for O2band = 1
        %Fn= figure(O2band+2) ;
        % set(Fn, 'Position', [488.2000 175.4000 739.2000 586.4000]);
        subplot(2+(O2band==1),2,k)
        hold on
        if k<3
            plot(Out.time(1:length(Out.FiFLD))*24,Out.FiFLD(:,O2band)*1E3,'o', 'MarkerEdgeColor',ColorWhite(1,:),'MarkerSize', MS); hold on
            plot(Out.time(10:length(Out.FiFLD)-10)*24,movmean(Out.FiFLD(10:end-10,O2band)*1E3,20),'m-','LineWidth',1);
        end
        plot(Out.time(1:end-1)*24,F(:,O2band)*1E3,'bx', 'MarkerSize', MS,'MarkerEdgeColor',[.5 .5 .5]); hold on
        plot(Out.time(10:end-11)*24,movmean(F(10:end-10,O2band),20)*1E3,'b-','LineWidth',1);
        


        plot(Out.time(2:length(Out.Fscope))*24,Out.Fscope(2:end,O2band),'Color', 'k');
        
        if k == 1
            legend('','iFLD','','BSF','SCOPE','location','south')
        end

        set(gca, 'ylim', [ -1 1.1], 'xlim', [6 18], 'xtick', 0:6:24)
        title([titles{k} O2bandname(O2band)])
        %xlabel('hour')
        if k ==1 || k == 3, ylabel('F (Wm^{-2}\mum^{-1}sr^{-1})'), end
        %     if d == 1, legend(['retrieved ' char(O2bandname(O2band))],['SCOPE '  char(O2bandname(O2band))] ), end
        if O2band==2 && k>2
            xlabel('hour')
            
        end
    end
end
%%
figure(6),clf
J = 494;

% Out = Noprior;
% for O2band = 1:2
%     switch O2band
%         case 1, data = Out(4).O2A(J);
%         case 2, data = Out(4).O2B(J);
%     end
%     
%     wl = data.wl;
%     E = data.E./data.normE;
%     piL = data.piL./data.normpiL;
%     piLr = data.piLr./data.normpiL;
%     a = Out(4).a(J);
%     
%     subplot(2,2,2*O2band-1)
%     plot(wl,[E, piL], 'x-', 'MarkerSize', 3);
%     xlabel('wl (nm)')
%     ylabel('norm. irradiance')
%     legend('E', 'piL', 'location', 'southeast')
%     
%     subplot(2,2,2*O2band)
%     hold on
%     plot(log(E), log(piL), 'kx','MarkerSize',3)
%     hold on
%     plot([-2.5 0],[-2.5 0], 'k--')
%     plot(log(E), log(piLr), 'r*', 'MarkerSize',3)
%     plot([-2.5 0],a*[-2.5 0],'r')
%     if O2band ==2, set(gca,'xlim', [-1 0],'ylim', [-1 0]),
%     else
%         set(gca,'xlim', [-3 0],'ylim', [-3 0]),
%     end
%     xlabel('log E(\lambda)/log(E_0)')
%     ylabel('log \piL(\lambda)/log(rE_0)')
% end
%%
Out = Noprior;
for O2band = 1:2
    
    figure(4+O2band),clf
    J = 494;
    for d = 1:2
        switch d
            case 2
                [~,J] = min(abs(Out(5).time-0.5));
                switch O2band
                    case 1
                        data = Out(5).O2A(J);
                    case 2
                        data = Out(5).O2B(J);
                end
            case 1
                [~,J] = min(abs(Out(10).time-0.5));
                %data = Out(10).O2A(J);
                switch O2band
                    case 1
                        data = Out(10).O2A(J);
                    case 2
                        data = Out(10).O2B(J);
                end
                
        end
        
        wl = data.wl;
        E = data.E./data.normE;
        piL = data.piL./data.normpiL;
        piLr = data.piLr./data.normpiL;
        a = data.a;
        
        spfig(2*d-1) = subplot(2,2,2*d-1); %#ok<*SAGROW>
        
        z=plot(wl,[E, piL], 'o-', 'MarkerSize', 4);
        set(z(2),'Marker','v')
        if d ==2
            xlabel('wl (nm)')
            
        end
        title(alphabet(2*d-1))
        ylabel('norm irradiance')
        legend('E', 'piL', 'location', 'southeast')
        switch O2band
            case 1, set(gca,'ylim',[0 1],'xlim',[758.8 770])
            case 2, set(gca,'ylim',[0 1],'xlim',[686.5 688.1])
        end
        
        spfig(2*d)=subplot(2,2,2*d);
        
        plot([-2.5 0],[-2.5 0], 'k')
        hold on
        plot(log(E), log(piL), 'ko','MarkerSize',4)
        plot(log(E), log(piLr), 'rv', 'MarkerSize',4)
        plot([-2.5 0],a*[-2.5 0],'r--','linewidth',1.5)
        % if d ==2, set(gca,'xlim', [-1 0],'ylim', [-1 0]),
        % else
        title(alphabet(2*d))
        switch O2band
            case 1,         set(gca,'xlim', [-2 0],'ylim', [-2 0]),
            case 2,         set(gca,'xlim', [-.6 0],'ylim', [-.6 0]),
        end
        % end
        if d ==2
            xlabel('log E(\lambda)/log(E_0)')
            
        end
        ylabel('log \piL(\lambda)/log(rE_0)')
    end
    
    resizefigure(spfig,2,2,.1,.13,.10,.08)
end