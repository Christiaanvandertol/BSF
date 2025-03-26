function MODTRAN = prepare_SRF_correction

%% loading SCOPE reflectances and fluorescence
p           = set_parameters;
SCOPEoutputfolder = 'OHP_2021-06-25-1444';
rsd         = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\rsd.csv'],',',2,0); %#ok<*DLMRD>
rdd         = dlmread(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\rdd.csv'],',',2,0);
wlS         = load(['..\output\SCOPE_simulation\', SCOPEoutputfolder, '\wlS.txt']);

%%
SZA         = [75,60,45,30];

%% load MODTRAN data
s1          = importdata('FLEX-S3_std.atm');%importdata('FLEX-S3_SPAIN_SZA75.atm');
s2          = importdata('FLEX-S3_SPAIN_SZA60.atm');
s3          = importdata('FLEX-S3_SPAIN_SZA45.atm');
s4          = importdata('FLEX-S3_SPAIN_SZA30.atm');
%%
for iSZA = 1:1%length(SZA)
    switch iSZA
        case 1, s = s1;
        case 2, s = s2;
        case 3, s = s3;
        case 4, s = s4;
    end

    % interpolations and pre-processing
    wl_MODTRAN  = transform_wvl_from_vac_to_air(s.data(:,2));
    index2      = find(wl_MODTRAN>p.wl_left(1) & wl_MODTRAN<p.wl_right(1));
    index3      = find(wl_MODTRAN>p.wl_left(2) & wl_MODTRAN<p.wl_right(2));

    T           = s.data(:,3:20);
    atmo.M      = [T(:,1) T(:,3) T(:,4) T(:,5) T(:,12) T(:,16)];
    atmo.Ta     = 25; % influence of air temperature is negligible, but a value
    % has to be provided.

    %interpolate SCOPE simulated reflectance to MODTRAN wavelengths
    SAIL.rsd    = interp1(wlS,rsd(5,:),wl_MODTRAN);
    SAIL.rdd    = interp1(wlS,rdd(5,:),wl_MODTRAN);

    % calculate irradiance
    [Esun,Esky] = calcIrradiance(atmo,SAIL,wl_MODTRAN);

    % interpolate SCOPE fluorescence to the MODTRAN wavelengths
    E_MODTRAN   = Esun;%+Esky;
    SAIL.rsd(isnan(SAIL.rsd)) = 0;

    % normalize the band depth by the interpolated values
    piL_MODTRAN0    = E_MODTRAN.*SAIL.rsd;

    normpiL2         = interp1([wl_MODTRAN(index2(1)),wl_MODTRAN(index2(end))],[piL_MODTRAN0(index2(1)),piL_MODTRAN0(index2(end))],wl_MODTRAN(index2));
    normpiL3         = interp1([wl_MODTRAN(index3(1)),wl_MODTRAN(index3(end))],[piL_MODTRAN0(index3(1)),piL_MODTRAN0(index3(end))],wl_MODTRAN(index3));

    MODTRAN.normpiO2A(:,iSZA)     = normpiL2;
    MODTRAN.normpiO2B(:,iSZA)     = normpiL3;
    MODTRAN.piL(:,iSZA)           = piL_MODTRAN0;
    MODTRAN.E(:,iSZA)             = E_MODTRAN;
    MODTRAN.Esky(:,iSZA)        = Esky;
    MODTRAN.Esun(:,iSZA)        = Esun;
    MODTRAN.wl(:,iSZA)         = wl_MODTRAN;
    MODTRAN.iO2A(:,iSZA)        = index2;
    MODTRAN.iO2B(:,iSZA)        = index3;
end
MODTRAN.SZA =  SZA;