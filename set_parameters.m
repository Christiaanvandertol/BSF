    function parameters = set_parameters

parameters.wl_in       = [761,686.9 715]; % next 6 coefficients are all used in iFLD
parameters.wl_out      = [755,686.5 721];
parameters.wl_left0    = [750,680 715];
parameters.wl_left1    = [755,686.5 715];
parameters.wl_right0   = [772,688.1 721];
parameters.wl_right1   = [777,690 735];
%parameters.wl_left     = [759, 686.5]; % just left, used in new method
parameters.wl_left     = [759, 686.1 715]; % just left, used in new method, 25 feb 2025
parameters.wl_right    = [768, 688.1 721]; % just right, used in new method

%parameters.wl_left     = [759, 685.0]; % left, used in new method, test 11 March 2025
%parameters.wl_right    = [768, 690.0]; % right, used in new method
