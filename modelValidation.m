clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ISP = 230;
g0 = 9.81;
g = 9.81;
mp = 72*5;
mf = 250;
m0 = mf+mp;
tb = 72;

ueq = ISP*g0;

vb = -ueq*log(mf/m0) - g*tb

hb = -ueq*tb*log(m0/mf)/((m0/mf) - 1) + ueq*tb - 0.5*g*tb^2

hmax = hb + vb^2/(2*g)