function [ wvl_a ] = transform_wvl_from_vac_to_air( wvl_v )

%
% :PURPOSE:
%
%   Function used to correct the wavelenghts of MODTRAN irradiance from vacuum to air, according to the formula of Edlén 1966
%
% :OUTPUTS:
% 
%   wvl_a: Corrected wavelenghts for the MODTRAN irradiance file
%    
%INPUT PARAMS:
%    
%    wvl_v: Wavelenghts of the different channels of the MODTRAN input irradiance file
%
%
%: Created:	25-feb-2013
%
% :AUTHOR:
% Sergio Cogliati
% sergio.cogliati@unimib.it



c=[0.99962627,5.4209431e-007,-1.2545200e-009,...
    1.4952310e-012,-8.9989458e-016,2.1672297e-019];

wvl_a = wvl_v .* ...
    (c(1) + c(2).*wvl_v + c(3).*wvl_v.^2 + c(4).*wvl_v.^3 + c(5).*wvl_v.^4 + c(6).*wvl_v.^5);



end

