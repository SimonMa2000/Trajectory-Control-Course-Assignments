function prob_04_studentName()
%
% This function computes a simulation of a the Lorenz system, creating
% a visualization of the Lorenz attractor fractal.
%

%~~~~~~~~~~~~~~~~~  Set up for the simulation  ~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Duration of 10 seconds, with a maximum time-step of 0.005 seconds.
% x(0) is randomly selected from [5, 35]
% y(0) is randomly selected from [-30, 5]
% z(0) is randomly selected from [-5, 35]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

clc;clear;
% Function randMy(a,b):
% generate a uniformly ditributed pseudo-random value in given range [a,b]
randMy = @(a,b)(a + (b-a).*rand());
% set up for the simulation
x0 = randMy(5, 35);
y0 = randMy(-30, 5);
z0 = randMy(-5, 35);
stateInit = [x0, y0, z0]';
tMax = 10;
hMax = 0.005;
tSpan = linspace(0,tMax,tMax/hMax+1);
dynFun=@(t,state)(LorenzDynamics(t,state));
% boundary check TBD
option = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);
[tOde45, zOde45] = ode45(dynFun, tSpan, stateInit, option);
tOde45 = tOde45'; zOde45 = zOde45';
%~~~~~~~~~~~~~~~~~~~~~~~~  Run the simulation  ~~~~~~~~~~~~~~~~~~~~~~~~~~~%
[time, state] = EulerMethodSimulation(dynFun, tSpan, stateInit, hMax);

%~~~~~~~~~~~~~~~~  Plot the trajectory in state space ~~~~~~~~~~~~~~~~~~~~%
% Use the plot3() command to plot (x(t) vs. y(t) vs. z(t)) on a single plot
% Set the axis commands so that the axis are equal in scale and not visible
% Add a title to the subplot
% The subplot() command places this plot on the right side of the figure.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Set up the figure
figure(1040); clf;
subplot(2, 3, [3,6]);hold on;
% plot the Lorenz attractor  (trajectory in state space)
plot3(state(1,:),state(2,:),state(3,:),'-r');
plot3(zOde45(1,:),zOde45(2,:),zOde45(3,:),'-b');
legend('Red - MyLorenzDynamics','Blue - ode45');
axis equal;
grid;
title('Lorenz attractor trajectory in state space x-y-z');
xlabel('X');
disp(['StateInit = [' num2str(stateInit') ']'])

%~~~~~~~~~~~~~~~  Plot the state as a function of time ~~~~~~~~~~~~~~~~~~~%
% The top sub-plot shows all three states as a function of time
% The bottom sub-plot shows the derivative of each state
% Use the legend command to label each curve
% Label both axis and add a title
% The subplot() commands split the left side of the figure into top and bottom.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h1 = subplot(2, 3, [1, 2]); hold on;

% plot the state vs time
plot(time,state(1,:),'-r');
plot(time,state(2,:),'-g');
plot(time,state(3,:),'-h');
xlabel('time (s)');
ylabel('state');
legend('x', 'y', 'z');
title('State vs Time');

% Plot the derivative of each state as a function of time:
h2 = subplot(2, 3, [4, 5]); hold on;
plot(time(2:end),diff(state(1,:)),'-r');
plot(time(2:end),diff(state(2,:)),'-g');
plot(time(2:end),diff(state(3,:)),'-h');
xlabel('time (s)');
ylabel('state');
legend('x', 'y', 'z');
title('Derivative of State vs Time');
%% FOR DEBUG: ode45 plot 
%{
% Plot the state as a function of time:
h1 = subplot(2, 3, [1, 2]); hold on;

%%%% TODO: plot the state vs time
plot(tOde45,zOde45(1,:),'-r');
plot(tOde45,zOde45(2,:),'-g');
plot(tOde45,zOde45(3,:),'-h');
xlabel('time (s)');
ylabel('state');
legend('x', 'y', 'z');
% Plot the derivative of each state as a function of time:

h2 = subplot(2, 3, [4, 5]); hold on;
%%%% TODO: plot the derivative of each state vs time
plot(tOde45(2:end),diff(zOde45(1,:)),'-r');
plot(tOde45(2:end),diff(zOde45(2,:)),'-g');
plot(tOde45(2:end),diff(zOde45(3,:)),'-h');
xlabel('time (s)');
ylabel('state');
legend('x', 'y', 'z');
%}
end
