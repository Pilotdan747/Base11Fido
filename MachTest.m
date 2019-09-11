clear
clc
close all


x = linspace(0,1);
x2 = linspace(1,4);

Cd = 0.33./sqrt(1-x.^2);
Cd2 = 0.33./sqrt(x2.^2-1);

figure
plot(x, Cd)
hold
plot(x2, Cd2)
