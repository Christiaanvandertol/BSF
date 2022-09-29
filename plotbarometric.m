function [hr,p_p02,costs] = plotbarometric

z = 0:100;
%p_p0 = barometric(z);

z1 = 100;

hr = 7:0.1:16;
ts = calczenithangle(265,hr,0,0,3+35/60,43+43/60);
costs = cos(ts);
costs0 = max(costs);
costo = 1;

p_p01 = barometric(z);%*(1+costs0));
%p_p02 = barometric(z1*(1+cos(ts/180*pi)));
p_p02 = barometric(z1);%*(1+costs));
%p_p03 = barometric(z1);%*(1+cos(50/180*pi)));

F4 = figure(4); hold on
set(F4,'Position',[458.6000 293 560 236.8000])

subplot(121)
%plot(z,1./(p_p01) + (1-p_p01)./p_p01.*costs0./costo);
%keyboard
%hold on
plot(z,1./p_p01.*(costs0/costo+1)-costs0/costo,'k');
%plot(z,1+(1./p_p01-1).*(costs0+1),'k');
%plot(z,1./p_p01.*(1+costs0)-costs0,'k');
xlabel('MH (m)')
ylabel('(x_0+x_1+x_2)/x_0')
%set(gca, 'ylim',[1 1.025])

subplot(122), hold on
plot(hr,1./p_p02 .* (costs/costo+1)-costs/costo,'k');
%plot(hr,1+(1./p_p02-1).*(costs+1),'k');
%hold on
%plot(hr,1./p_p03 .* (
xlabel('time (hour of day)')
ylabel('(x_0+x_1+x_2)/x_0')
set(gca, 'ylim',[1 1.025])
