%% IWP MQP Full Scale Thrust-Pulsing Optimization
clc; close all; clear all

% Airship parameters
rho = 1.225; 
muair = 1.79E-5;

L = 11.0983; 
d = 2.3121;
m = 266.6999/9.81;
vol23 = 9.8823;
vol23ft = 106.3718204;
FR = 4.8;
FF_body = 1.2059;
FF_tail = 1.2306; 
SA_body = 64.19736;
SA_tail = 3.146721;
C_tail_avg = 0.844798; 

% Operational conditions
% wind direction input
disp('Wind Direction (in degrees):');
disp(' 0  --> Wind pushing the airship forward (Wind from behind)');
disp(' 90  --> Wind from the right (starboard side)');
disp(' 180 --> Wind blowing directly into the airship (headwind)');
disp(' 270 --> Wind from the left (port side)');

wind_angle = input('enter the wind direction in degrees (0 to 360): ');

while wind_angle < 0 || wind_angle > 360
    disp('Invalid input! Wind direction must be between 0 and 360 degrees.');
    wind_angle = input('enter the wind direction in degrees (0 to 360): ');
end

vw = 5; % average airspeed

v_0 = 4; % Initial velocity at t=0 in m/s

vinf = @(v) v - vw*cosd(wind_angle);

K = 1.3118; % drag due to lift factor (K)

% dynamic pressure
q = @(v) 0.5*rho*(vinf(v)^2); 

% Skin friction drag coefficients
Cf_body = 0.003661479;
Cf_tail = 0.005968424;

% Drag coefficients
CD0_body = 0.028684018;
CD0_tail = 0.002338768;

CD0_gond = 0.075485467;
CD0_eng = 0.039954191;
CD0_engmount = 0.009911004;
CD0_cable = 0.09617812;
CD0_gear = 0.008667059;
CD0_int = 1.50264E-05;

CD0 = CD0_body + CD0_tail + CD0_gond + CD0_eng + CD0_engmount + CD0_cable + CD0_gear + CD0_int;

L_aero = 20.00253758;
CL_aero = 0.206539206;

% Total drag
F_drag = @(v) (CD0 + K*(CL_aero^2))*q(v)*vol23*sign(vinf(v));

% IWP thrust function
F_thrust = @(v) 21.67642589 * (1 - heaviside(v - 3.36));

% Define the differential equation
odefun = @(t, v) (F_thrust(v) - F_drag(v)) / m;

% Time span and initial condition for ode45
tspan = [0 900]; % Time from 0 to 900s (15m)
initial_conditions = v_0; % Initial velocity at t=0

% Solve ODE using ode45
[t, v] = ode45(odefun, tspan, initial_conditions);

% Calculate F_drag and F_thrust for the entire time interval
F_drag_values = arrayfun(@(v) F_drag(v), v);
F_thrust_values = arrayfun(@(v) F_thrust(v), v);
F_q_values = arrayfun(@(v) q(v), v);

% Plot results
figure(1);
plot(t, v, '-','LineWidth',1.2);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
ylim([-5, 5]);
xlim([0, 60]);
title('Airship Velocity over Time');
grid on;

% Plot F_drag
figure(2);
plot(t, F_drag_values, '-r', 'LineWidth', 1.2);
hold on
plot(t, F_thrust_values, '-b', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Force (N)');
xlim([0, 60]);
title('Drag vs. IWP thrust over Time');
grid on;
legend('Drag', 'IWP Thrust')

