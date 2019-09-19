% Trajectory Analysis
% Author: Daniel Owen
% Created: 1/25/19
%
% Version 1.0 - 2/2/19
%   Working Simulation - See TODO for improvements
%
% From MIT PDF:
%   http://web.mit.edu/16.unified/www/FALL/systems/Lab_Notes/traj.pdf
%
% t - Time                        Cd - Drag Coefficient
% h - altitude                    A - reference area
% V - velocity pos up             mdotf - fuel mass flow
% F - force pos up                ue - exhuast velocity
% D - drag                        i  - time index
% T - propulsive trust
% dt - time step
% rho - density
% g - gravity
% m - mass
%
% hi+1 = hi + Vi*(ti+1 - ti)
% Vi+1 = Vi + (-g - (1/2)*rho*Vi*abs(Vi)*(Cd*A/mi) +
% (Vi/abs(Vi))*(mdotfi*uei/mi)) * (ti+1 - ti)
% mi+1 = mi + (-mdotfi) * (ti+1 - ti)
%
%
% Assumptions
%   Cd = 0.550 (Saturn V)           NEEDS FIXING
%   ISP = 309.78 s
%   A = 0.25^2*pi (0.5m diameter)
%   mdot = 2.75 kg/s
%   burn time = 106.49 s
%   mass = 500 kg
%
% Airbrakes
%   A increases by 4 time
%   Cd increases to 1.2 (flat plate)
%
% TODO
%   Add Rho as a function of height - Done
%   Add Cd as a function of V (maybe as a function of Mach)
%   Model to get more acuurate Cd - Gage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INIT
%Constants
g = 9.81; %m/s^2 

%Atmosphere model
[Z, Z_L, Z_U, T, P, rho, c, g1, mu, nu, k, n, n_sum] = atmo(200, 0.01, 1);
rho1 = rho';
P1 = P';
c1 = c';

%Assumptions
Cd = 0.33;
dBody = 0.35;
ILimit = 889600; 

% Engine Params
ISP = 230;
mdot = 5.0;
dryMass = 250;       %CHECK THIS
dExhaust = 0.14783;
exitArea = dExhaust^2/4*pi;   %m^2
ue = ISP*g;     
pExhaust = 1e5;

%IFactor = 1;

A = dBody^2/4*pi;              %Exhaust Velocity
propMass = ILimit/ue;
mass = propMass + dryMass;
burnTime = propMass/mdot;

dV = ue*log(mass/dryMass);

%%
% thrust = mass*g*1.73;
% mdot = thrust/ue;
% burnTime = propMass/mdot;
%%

%Sim time in seconds - Arbitrally large - cut down data size later
maxTime = 1000000;

%Find array size based on time increment
dt = 0.25; %seconds
iMax = floor(maxTime/dt);

%State Varaibles
mdotf = zeros(1,iMax+1);
h = zeros(1,iMax+1);       %Height
V = zeros(1,iMax+1);       %Velocity
m = zeros(1,iMax+1);       %Mass
Cdm = zeros(1,iMax+1);     %Drag Coefficient
thrust = zeros(1,iMax+1);  %Thrust
q = zeros(1,iMax+1);       %Dynamic Pressure

drag = zeros(1,iMax+1);    %Drag force
M = zeros(1,iMax+1);       %Mach number

%Impulse Check
I = zeros(1,iMax+1);

%Inital Conditions
h(1) = 1371;         %4500 feet/ 1371 meters
V(1) = 0.0001;       %nonzero small number to avoid devide by 0 errors
m(1) = mass;

%Define mdot for burn time and 0 all else
for i = 1:(burnTime/dt)
    mdotf(i) = mdot;
end

%% INTIGRATION
%Integration Initilization
i = 1;
count = 0;
flag = 1;

%Integration - Loop until flag set to 0
while flag ~= 0
    %Atmosphere at time step
    atmIndex = max(floor(h(i)/10), 1);
    
    %Density, pressure, and speed of sound
    rho = rho1(atmIndex);
    p = P1(atmIndex);
    if atmIndex > 86/0.01
        a = inf;
    else
        a = c1(atmIndex);
    end
    
    %Aerodynamic state at timestep
    q(i) = (1/2)*rho*V(i)*abs(V(i));
    M(i) = V(i)/a;
    
    %Forces at time step
    
    %min(3, stuff) account for quicker transonic crossing
    %Sets max Cdm to 3 - used to blow up to inifinty
    if (abs(M(i)) < 1)
        Cdm(i) = min(3,Cd/sqrt(1 - M(i)^2));
    else
        Cdm(i) = min(3,Cd/sqrt(M(i)^2 - 1));
    end
    
    drag(i) = q(i)*(Cdm(i)*A/m(i));  %Drag accel (Not really a force)
    
    if (mdotf(i) ~= 0)
        thrust(i) = V(i)/abs(V(i))*((mdotf(i)*ue + (pExhaust - p)*exitArea)/m(i));
    else
        thrust(i) = 0;
    end
    
    %MIT Equations
    h(i+1) = h(i) + V(i)*dt;
    
    V(i+1) = V(i) + (-g-drag(i)+thrust(i)) * dt;
    
    m(i+1) = m(i) + -mdotf(i)*dt;
    
    
    %Impulse Integrator
    I(i+1) = I(i) + thrust(i) * m(i) * dt;
    
    % Deploy airbrakes at Apogee
    if (h(i+1) < h(i) )   %&& h(i) < 10000)
        A = 4*0.25^2*pi;
        Cd = 1.28;
        
        %End integration when rocket lands at 4500 feet/ 1371 meters
        if (h(i) < 1371) 
            flag = 0;
        end
    end
    
    i = i + 1;
end

%% PLOTING
%Reduce Data Size down to flight time (liftoff to landing)
h = h(1:i);
V = V(1:i);
m = m(1:i);
drag = drag(1:i);
M = M(1:i);
Cdm = Cdm(1:i);
I = I(1:i);
thrust = thrust(1:i);
q = q(1:i);

%Display total impulse
I(i)

%Plotting Time
t = 0:dt:(i*dt)-dt;

%Hight by time
figure
plot(t, h)
xlabel('Time (s)')
ylabel('Height (m)')
title('Height Plot')
line([0,t(length(t))],[1e5,1e5])

%Velocity by time
figure
plot(t, V)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Plot')

%Drag by time
figure
plot(t, drag.*m)
xlabel('Time (s)')
ylabel('Drag (N)')
title('Darg Plot')

%Mach Number by time
figure
plot(t, M)
xlabel('Time (s)')
ylabel('Mach Number')
title('Mach Plot')

%Accel Plot
figure
plot(t, (thrust - drag - g)./g)
xlabel('Time (s)')
ylabel('Acceleration (Gs)')
title('Acceleration Plot')

%Thrust by time
figure
plot(t, thrust.*m/1000)
xlabel('Time (s)')
ylabel('Thrust (kN)')
title('Thrust Plot')

%Thrust by time
figure
plot(t, q)
xlabel('Time (s)')
ylabel('q (Pa)')
title('Dynamic Pressure Plot')
