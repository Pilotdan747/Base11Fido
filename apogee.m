% Apogee Equation
% Author: Daniel Owen
% Created: 1/24/19
%
% Impulse Limit - 889,600

mass = 500;

thrust = [5e3:1:30e3];
impulse = 889600;
time = impulse./thrust;

x = 1/2*(thrust./mass - 9.81).*time.^2;

figure
plot(thrust, x);
xlabel('Thrust (N)')
ylabel('Altitude at Burn Out (m)')
title('Burn Out Altitude')

v = (thrust./mass - 9.81).*time;

figure
plot(thrust, v);
xlabel('Thrust (N)')
ylabel('Velocity at Burn Out (m/s)')
title('Burn out velocity')

t2 = v./9.81;

x2 = x + v.*t2 + 0.5*-9.81.*t2.^2;

figure
plot(thrust, x2);
xlabel('Thrust (N)')
ylabel('Apogee (m)')
title('Apogee')
