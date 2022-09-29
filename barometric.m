function p_p0 = barometric(z)

g = 9.81;
M = 0.02896968; % kg/mol
T0 = 310;%288.16; %K
cp = 1004; %J/(kg K);
R0 = 8.314462618; %J/(mol K)

p_p0 = (1-g*z/(cp*T0)).^(cp*M/R0);