

k = 4;
O2bandname ={'O2A','O2B'};
for O2band = 1:2
switch O2band
    case 1,     O2 = Out(4).O2A; limit = 57;
    case 2, O2 = Out(4).O2B; limit = 9;
end

[Ri,Ji]= deal(zeros( length(O2(1).RES), length(O2)));
[sF,F] = deal(zeros(length(O2),1));
for j = 1:length(O2)
    R = O2(j).RES(1:limit);
    J = O2(j).JAC(1:limit);
    Ri(:,j) = R;
    Ji(:,j) = J;
    sF(j) = sqrt(inv(J'*J)*var(R));
    F(j) = O2(j).F;
end
sF = sF*1E3;
F = F*1E3;


%%
figure(1), 
% subplot(131)
% z = errorbar( Out(4).time(1:length(F))*24,F,sF);
% % z = plot(Out(4).time(1:length(F))*24,[F-sF,F, F+sF],'ko','MarkerSize',3);
% % set(z(1), 'MarkerEdgecolor',[.7 .7 .7])
% % set(z(3), 'MarkerEdgecolor',[.7 .7 .7])
% ylabel('F +/- \sigma F (W m^{-2}\mum^{-1}sr^{-1})')
% xlabel('time (hour)')

subplot(2,2,O2band*2-1)
z = plot( Out(4).time(1:length(F))*24,F);
ylabel('F (W m^{-2}\mum^{-1}sr^{-1})')
xlabel('time (hour)')
title(O2bandname{O2band})

subplot(2,2,O2band*2)
z = plot( Out(4).time(1:length(F))*24,sF);
ylabel('\sigmaF (W m^{-2}\mum^{-1}sr^{-1})')
xlabel('time (hour)')
title(['uncertainty ' O2bandname{O2band}])

end