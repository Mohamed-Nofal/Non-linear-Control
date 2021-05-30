clc;clear all
zeta = 1.1547;
w_n  = 3.4641;
ts = 0.5;

S_1 = -zeta*w_n+i*w_n*sqrt(1-zeta^2)
S_2 = -zeta*w_n-i*w_n*sqrt(1-zeta^2)

Z_1 = exp(ts*S_1)
Z_2 = exp(ts*S_2)

z = tf('z',ts)
Am= tf((z-Z_1)*(z-Z_2),ts,'Variable','z^-1')