function [E,a,y2,piLmod] = cost4F(F,input)
%F=0;%max(0,F);
logx        = input.logx;
y           = input.y;
cos_sza     = input.cos_sza;
cos_vza     = input.cos_vza;
fwlf        = input.flwf;
normpiL     = input.normpiL;% 
SRC         = input.SRC;
%residual    = input.residual;
%y(y==0)     = 0.01; % In principle this line is not necessary but 
%                       in my 1-yr time series of data, there were a few 
%                       occasions where at one wl, the radiance was zero (recording error)
%                       to prevent the code from stopping, I added this
%                       line.
priorweight = input.priorweight;

if input.atcor
    a = input.aprior;
    I = find(logx<input.logxlim); % only use the deeper part of the band for fitting
    % a loop is needed due to the reabsorption of F by O2 on the path from the
    % surface to the sensor. This interferes with the estimate of the path
    % length a. Three steps is sufficient.
    for k = 1:3
        % The fluorescence is subtracted from the normalized radiance to
        % obtain a pure reflected radiance. The fluorescence varies
        % linearly. It must also be normalized by normpiL, and reabsorption
        % is calculated.        
        
        Fra = F*fwlf.* exp(logx*(a-1)./(1+cos_vza/cos_sza));
        keyboard
        y2      = (y.*normpiL - Fra)./(normpiL - F*fwlf);%
%        y2      = (y.*normpiL - Fra)./(normpiL - F);%
        
        y2(y2<0) = 1E-12;  
        logy2   = log(y2) -SRC(1:length(y2))*( (a-1));%- residual;

        if priorweight>0
            a       = input.aprior;%logx \ logy2; % this is linear regression
        else
            a       = logx \ logy2; % this is linear regression
        end
    end
else
    a = input.a;
    %I = find(logx<input.logxlim); % use the whole band for fitting.  
    [~,I] = min(logx);
    % Because the path length 'a' is already known, it is not necessary to iterate here.
    Fra = F*fwlf.* exp(logx*(a-1)./(1+cos_vza/cos_sza));
    y2      = (y.*normpiL - Fra)./(normpiL - F*fwlf);%
    logy2   = log(y2) - SRC*( (a-1));
end
logymod     = a*logx; % this is the modelled radiance
E           = 1E3*[logymod(I)-logy2(I)+ priorweight*(a - input.aprior)]; % this is the cost function
E           = E(~isnan(E));

piLmod       = y2/a.*normpiL+F*fwlf;