w_F     = ones(size(Out.piL0fluo,1),1);

w_F(Out.wlfluo>690 & Out.wlfluo<735 | Out.wlfluo>770) = 0.2;

%w_F(Out.wlfluo>735 & Out.wlfluo<755) = 1E3;



opt_alg   = 'tr';    % optimization algorithm
stio      = 'off';   % options:


J = find(Out.wlfluo>650 & Out.wlfluo<800);
clear Fall
for k = 1:size(Out.E0fluo,2)
    [x_F,f_wvl_F,r_wvl_F,resnorm_F,exitflag_F,output_F]  = FLOX_SpecFit_6C(Out.wlfluo(J),Out.E0fluo(J,k),Out.piL0fluo(J,k),[1,1],w_F(J),opt_alg,stio,Out.wlfluo(J));
    Fall(:,k) = f_wvl_F(:);
    rall(:,k) = r_wvl_F(:);
end
