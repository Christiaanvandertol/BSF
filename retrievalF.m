function [O2A,O2B] = retrievalF(wl,E,piL,opt,aprior,cos_sza,cos_vza,priorweight,p,SRCA,SRCB,flwf,residual)
 
%% loop over the bands
for O2band = 1:2
    % select the band and interpolate the (ir)radiances over the band.
    index2          = find(wl>p.wl_left(O2band) & wl<p.wl_right(O2band));

    switch O2band
        case 1 
            I = find( (wl>755& wl<759) | (wl>770 & wl<775));
        case 2
            I = [index2(1)-20:index2(1),index2(end):index2(end)+20];
            
    end

    %different ways to do the following interpolation. The first one works
    %best for synhtetic and real data.
    normE       	= interp1(wl(I),E(I),wl(index2));
    x = polyfit(wl(I),piL(I)./E(I),2);
    normr = polyval(x,wl(index2));
    normpiL         = normE.*normr;
   
  %  normE       	= interp1([wl(index2(1)),wl(index2(end))],[E(index2(1)),E(index2(end))],wl(index2));
  %  normpiL         = interp1([wl(index2(1)),wl(index2(end))],[piL(index2(1)),piL(index2(end))],wl(index2));

   % normE       	= interp1([wl(index2(1)-20:index2(1));wl(index2(end):index2(end)+20)],[E(index2(1)-20:index2(1));E(index2(end):index2(end)+20)],wl(index2), 'spline');
   % normpiL       	= interp1([wl(index2(1)-20:index2(1));wl(index2(end):index2(end)+20)],[piL(index2(1)-20:index2(1));piL(index2(end):index2(end)+20)],wl(index2), 'spline');

   % normE       	= interp1([wl(index2(I1)),wl(index2(end))],[E(index2(I1)),E(index2(end))],wl(index2), 'linear', 'extrap');
  %  normpiL         = interp1([wl(index2(I1)),wl(index2(end))],[piL(index2(I1)),piL(index2(end))],wl(index2), 'linear', 'extrap');

  % normE       	= interp1([wl(index2(1)),wl(index2(end))],[mean(E(index2(1)-10:index2(1))),mean(E(index2(end):index2(end)+10))],wl(index2));
   % normpiL         = interp1([wl(index2(1)),wl(index2(end))],[mean(piL(index2(1)-10:index2(1))),mean(piL(index2(end):index2(end)+10))],wl(index2));

    %normE       	= ones(length(wl(index2)),1) * E(index2(1));
    %normpiL         = ones(length(wl(index2)),1) * piL(index2(1));


    % normalize the band depth by the interpolated values
    input.logx      = log(E(index2)./normE);
    input.y         = piL(index2)./normpiL;

    input.normpiL 	= normpiL;
    input.aprior    = aprior;
    input.cos_sza   = cos_sza;
    input.cos_vza   = cos_vza;
    input.priorweight = priorweight;
    
    % specific settings per band (O2A and O2B)
    switch O2band
        case 1 %O2A

            if nargin>11
                input.flwf = flwf.O2A;
            else
                input.flwf      = 0.7+0.3*(length(index2):-1:1)'/length(index2);
            end
            % fluorscence varies approximately linearly over the band, with
            % a 34 percent reduction between 759 and 768 nm.
            %input.flwf      = 0.66+0.34*(length(index2):-1:1)'/length(index2);
            %input.flwf      = 0.7+0.3*(length(index2):-1:1)'/length(index2);
            input.atcor     = 1;
            input.logxlim   = 0;% -0.7 
            input.SRC       = SRCA;
        case 2 %O2B
            % I do the fitting of the opticalp path lenght only for the O2A band, not the
            % O2B. For the O2B band I adopt the fit obtained with O2A.
            % This is because the O2B filling is much less strong, it is
            % easier to use the O2A for this purpose.
            %input.flwf      = 1;
            input.flwf      = 1-0.1*(length(index2):-1:1)'/length(index2);
            input.atcor     = 0;
           % input.atcor     = 1;
           input.SRC = SRCB(1:length(index2));
            if ~isnan(a)    
                input.a = a;
            else
                input.a = 1;
            end
            input.logxlim = 0;%-0.3;
 %           keyboard
    end
    
    % the following iteratively changes fluorescence until the depth of the
    % irradiance and reflected radiance is linearly related.
    f           = @(Fi)cost4F(Fi,input);
    if min(input.logx)<input.logxlim 
        %[F,RESNORM,RESIDUAL,EXITFLAG]           = lsqnonlin(f,0,-1E10,1E10,opt); %#ok<ASGLU> % F normalized by the radiance outside the band
        [F,RESNORM,RESIDUAL,EXITFLAG]           = lsqnonlin(f,0,-1E10,1E10,opt); %#ok<ASGLU> % F normalized by the radiance outside the band
        
       % if O2band ==1 , test_sensitivity(input), end
        
    else
        F = 0; EXITFLAG = 1;
    end
    
    [~,a,piLr,piLmodb]       = cost4F(F,input); % the slope of the regression, which is the ratio of optical depths
%keyboard
%    piLmod = piL;
 %   piLmod(index2)= piLmodb;
  %  [F] =  iFLD(E,piLmod,wl,p.wl_left0(O2band),p.wl_left1(O2band),p.wl_in(O2band),p.wl_right0(O2band),p.wl_right1(O2band),p.wl_out(O2band));


    %% iFLD 
    iFLDr       =  iFLD(E,piL,wl,p.wl_left0(O2band),p.wl_left1(O2band),p.wl_in(O2band),p.wl_right0(O2band),p.wl_right1(O2band),p.wl_out(O2band));
    
    %% assign the output to the output structures
    O2.wl = wl(index2);
    O2.F = F;% fluorescence, the mean over the band (improve if more accuracy is needed)
    O2.E = exp(input.logx).*normE;%E(index2)./normE;
    O2.piL = input.y.*normpiL;
    O2.piLr = piLr.*normpiL;
    O2.normpiL = normpiL;
    O2.normE = normE;
    O2.a = a;
    O2.EXITFLAG =  EXITFLAG;
    O2.iFLD = iFLDr;
    O2.RESIDUAL = RESIDUAL;
   % keyboard

    switch O2band
        case 1
            O2A = O2;
        case 2     
            O2B = O2;
    end
end