%%  FLD and 3FLD and iFLD
%  SIF-Box data Process
%  Created by PY(p.yang@utwente.nl)
%  Date: 2016-10-05
%  E is irradiance/pi, L is radiance, wl is wavelength vector, wl_in is the
%  wavelength of inband (for example 760)

%% reference: 
% Alonso, L., Gómez-Chova, L., Vila-Francés, J., Amorós-López, J., Guanter, L., Calpe, J., & Moreno, J. (2008). 
% Improved Fraunhofer Line Discrimination method for vegetation fluorescence quantification. IEEE Geoscience and Remote Sensing Letters, 5(4), 620-624.
function [F] =  iFLD(E,L,wl,wl_left0,wl_left1,wl_in,wl_right0,wl_right1,wl_out0)
% find the in, left and right band location in the wl file
nwl                 =   length(wl);
if size(E,1)    ~= nwl
    E =  E';
end 
if size(L,1)    ~= nwl
    L =  L';
end 
[~,index_left0]     =   min(abs(wl-wl_left0));
[~,index_left1]     =   min(abs(wl-wl_left1));
[~,index_in]        =   min(abs(wl-wl_in));
[~,index_out0]       =   min(abs(wl-wl_out0));

[~,index_right0]    =   min(abs(wl-wl_right0));
[~,index_right1]    =   min(abs(wl-wl_right1));

index_left          =   index_left0:index_left1;
index_right         =   index_right0:index_right1;
index_out           =   [index_left0:index_left1,index_right0:index_right1];

wl_left             =   wl(index_left);
wl_right            =   wl(index_right);
wl_in               =   wl(index_in);
wl_out              =   [wl_left;wl_right];
% keyboard
% [L_out,index_left_i       =   max(L(index_left,:));

Eout_range          =   E(index_out,:);
Eout                =   E(index_out0,:);
Ein                 =   E(index_in,:);
Lout                =   L(index_out0,:);
Lin                 =   L(index_in,:);

Refl                =   L./E;

aR                  =   Refl;                                       % apparent reflectance
tR_out_range        =   aR(index_out,:);                            % true reflectance outside of the absorption line
tR_in               =   interp1(wl_out,tR_out_range,wl_in);               % tune reflectance inside the absoprtion line
E_in_interp         =   interp1(wl_out,Eout_range,wl_in);             % tune reflectance inside the absoprtion line


coeff_R             =   aR(index_out0,:)./tR_in;
coeff_F             =   Eout./E_in_interp.*coeff_R;
F                   =  ((coeff_R.*Eout).*Lin-Ein.*Lout)./((coeff_R.*Eout)-(coeff_F.*Ein));