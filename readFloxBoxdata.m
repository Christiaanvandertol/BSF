%clear all,%
close all
%files = dir('L*.csv');

%datapath = '180914';
%datapath = '180923';
%datapath = 'OHP/180926';
datapath = 'OHP/180927';
%datapath = 'IT1/180713';
%datapath = 'SWI/180804';

colorsFLD = {'r','m'};
colorsNEW = { 'k','b'};

%% the following indices are needed for the iFLD method
parameters  = set_parameters;

wl_in       = [761,686.9]; % next 6 coefficients are all used in iFLD
wl_out      = [755,686.5];
wl_left0    = [750,680];
wl_left1    = [755,686.5];
wl_right0   = [772,688.1];
wl_right1   = [777,690];
wl_left     = [759, 686.5]; % just left, used in new method
wl_right    = [768, 688.1]; % just right, used in new method


%%
for FB = 1:2%2
    switch FB
        case 1
            files = dir([datapath '/*_*09*.csv']);
            i_E = 1;
            i_piL = 5;
            %First = 1;
            Last = 210;
        case 2
            files = dir([datapath '/L*.csv']);
            i_E = 5;
            i_piL = 7;
            %First = 60;
            %Last = 700;
            Last = 600;
            %Last = 292;
    end
    
    for k = 1:length(files)
        filename = [datapath '/' files(k).name];
        fileID = fopen(filename);
        T = readtable(filename, 'TreatAsEmpty', '#N/D');
        fclose(fileID);
        x = table2array(T);
        D(k).filename = filename;
        D(k).wl = x(:,1);
        D(k).data = x(:,2:Last+1);
        
        format = '%4c';
        for n = 1:size(x,2)
            format =[format '%2c%d%c%d%c%d%c']; %#ok<AGROW>
        end
        
        fid = fopen(filename);
        line = fgetl(fid);
        z = sscanf(line,format);
        fclose('all');
        hr = z(7:8:end);
        minute = z(9:8:end);
        sec = z(11:8:end);
        D(k).time = datenum(0,0,0,hr,minute,sec); %#ok<*SAGROW>
    end
    
    %%
    %close all
    [a] = deal(nan*zeros(Last,1));
    [F,FY,iFLDr] = deal(nan*zeros(Last,2));
    
    %for I = 1:Last%length(F)%First : Last
     for I = 1:length(F)%First : Last
        wl          = D(i_E).wl;
        E           = D(i_E).data(:,I);
        piL         = D(i_piL).data(:,I);
        
       
        
        for O2band = 1:2
            index2  = find(wl>wl_left(O2band) & wl<wl_right(O2band));
            normE = interp1([wl(index2(1)),wl(index2(end))],[E(index2(1)),E(index2(end))],wl(index2));
            normpiL = interp1([wl(index2(1)),wl(index2(end))],[piL(index2(1)),piL(index2(end))],wl(index2));
            
            input.logx         = log(E(index2(1:end))./normE);
            input.y         = piL(index2(1:end))./normpiL;
            input.normpiL   = normpiL;
            switch O2band
                case 1 %O2A
                    input.flwf      = 0.6+0.4*(length(index2):-1:1)'/length(index2);
                    input.atcor = 1;
                case 2 %O2B
                    input.flwf = 1;
                    input.atcor = 0;
                    input.a = ai;
            end
            
            f               = @(Fi)cost4F(Fi,input);
            %FY(I,O2band)       =    lsqnonlin(f,0,0,.1);
          
            
            %F(I,O2band)       =    lsqnonlin(f,0,0,.1);
            F(I,O2band)       =    lsqnonlin(f,0,0,20);
            
            %[e,ai]        = cost4F(FY(I,O2band),input);

            [e,ai]        = cost4F(F(I,O2band),input);
            %F(I,O2band)        = mean(normpiL.*input.flwf)* FY(I,O2band);
            piLr = piL - F(I,O2band);
            
            a(I)        = ai;
            %F(I,O2band)        = piL(index2(1))* FY(I,O2band)*mean(input.flwf);
            

            %if O2band ==1 % I only apply iFLD to the O2A band
            [iFLDr(I,O2band)] =  iFLD(E,piL,wl,wl_left0(O2band),wl_left1(O2band),wl_in(O2band),wl_right0(O2band),wl_right1(O2band),wl_out(O2band));
            %iFLDr(I,O2band) = FLD3(E,piL,wl,wl_left(O2band),wl_in(O2band),wl_right(O2band));
            %iFLDr(I,O2band) = FLD(E,piL,wl,wl_left(O2band),wl_in(O2band));
            %end
            
            %  toc
            if I == 100
                figure(O2band),
                subplot(2,2,2*FB-1)
                plot(wl(index2),[E(index2)./normE piL(index2)./normpiL], 'x-', 'MarkerSize', 3);
                xlabel('wl (nm)')
                ylabel('norm. irradiance')
                legend('E', 'piL', 'location', 'best')
                %   keyboard
                
                
                subplot(2,2,2*FB)
                hold on
                plot(log(E(index2)./normE), log(piL(index2)./normpiL), 'kx','MarkerSize',3)
                hold on
                plot([-2.5 0],[-2.5 0], 'k--')
                %   if O2band ==1
                plot(log(E(index2)./normE), log(piLr(index2)./normpiL), 'r*', 'MarkerSize',3)
                %  end
                plot([-2.5 0],a(I)*[-2.5 0],'r')
                if O2band ==2, set(gca,'xlim', [-1 0],'ylim', [-1 0]), else
                    set(gca,'xlim', [-3 0],'ylim', [-3 0]), end
            end
            
        end
    end
    
    
    
    xlabel('log E(\lambda)/log(E(759))')
    ylabel('log \piL(\lambda)/log(\piL(759))')
    %%
    figure(3)
    for O2band = 1:2
        %figure(O2band+2),
        
        subplot(2,3,3*FB-2)
        plot(D(k).time(1:Last)*24,a,'k')
        ylabel('(x2+x3)/x1'),
        hold on
        plot(D(k).time(1:Last)*24,a*0+1,'k')
        %set(gca,'ylim', [0.98, 1.08])
        set(gca, 'xlim', [7,17])
        xlabel('hour')
        
        subplot(2,3,3*FB-2+O2band),
        plot(D(k).time(1:Last)*24,iFLDr(:,O2band)*1E3,'r'); hold on
        plot(D(k).time(1:Last)*24,F(:,O2band)*1E3,'k');
        ylabel('F (W m^{-2}\mu m^{-1}sr^{-1})');
        switch O2band
            case 1, set(gca, 'ylim', [-0.5 1.0])
            case 2, set(gca, 'ylim', [-0.50 1.0])
        end
        set(gca, 'xlim', [7,17])
        legend('F_{iFLD}','F_{new}', 'Location','North','orientation', 'horizontal')
        xlabel('hour')
        
        %    subplot(2,3,3*FB),
        %   plot(D(k).time(1:Last)*24,FY(:,O2band)*1E3,'k');
        %   ylabel('F/E');
        %   %set(gca, 'xlim', [D(k).time(First)*24 D(k).time(Last)*24])
        %   set(gca, 'xlim', [7,17])
        %   % set(gca, 'ylim', [0 1.2])
        %  xlabel('hour')
    end
end
% end

%%

