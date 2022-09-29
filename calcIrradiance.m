function [Esun, Esky] = calcIrradiance(atmo,SAIL,wl)

% Extract MODTRAN atmosphere parameters at the SCOPE wavelengths
 t1  = atmo.M(:,1);
 t3  = atmo.M(:,2);
 t4  = atmo.M(:,3);
 t5  = atmo.M(:,4);
 t12 = atmo.M(:,5);
 t16 = atmo.M(:,6);

% radiation fluxes, downward and upward (these all have dimenstion [nwl]
% first calculate hemispherical reflectances rsd and rdd according to SAIL
% these are assumed for the reflectance of the surroundings
% rdo is computed with SAIL as well
rsd = SAIL.rsd;
rdd = SAIL.rdd;
% assume Fd of surroundings = 0 for the momemnt
% initial guess of temperature of surroundings from Ta;

Fd      = zeros(length(rsd),1);
Ls      = equations.Planck(wl,atmo.Ta+273.15);

% Solar and sky irradiance using 6 atmosperic functions
%keyboard
Esun   = pi*t1.*t4;
Esky   = pi./(1-t3.*rdd).*(t1.*(t5+t12.*rsd)+Fd+(1-rdd).*Ls.*t3+t16);