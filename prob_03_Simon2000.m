function prob_03_SaiMa()
%
% This function computes a simulation of a simple pendulum.
%

%~~~~~~~~~~~~~~~~~  Set up for the simulation  ~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Duration of 20 seconds, with a maximum time-step of 0.01 seconds.
% Initial angle is randomly selected from [-3, 3] radians
% Initial angular rate is randomly selected from [-1, 1] radians / second
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%% Set up for the simulation
clc;clear;

param.g = 9.8;
param.l = 2;
nGrid = 200;
tGrid= linspace(0, 8, nGrid);
z0 = [0.01; 0]
dynFun=@(t,z)(PendulumDynamics(t,z));
hMax=tGrid(2)-tGrid(1)
%%  Simulation with ode45 and my function
option = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
%~~~~~~~~~~~~~~~~~~~~~~~~  Run the simulation  ~~~~~~~~~~~~~~~~~~~~~~~~~~~%

[tOde45, zOde45] = ode45(dynFun, tGrid, z0, option);
tOde45 = tOde45'; zOde45 = zOde45';

[tGrid, zEuler] = EulerMethodSimulation(dynFun, tGrid, z0, hMax)

%~~~~~~~~~~~~~~~~~~~  Make plots of the simulation  ~~~~~~~~~~~~~~~~~~~~~~%
% Create a single figure with three sub-plots (three rows, one column)
% The top sub-plot is pendulum angle vs time
% The middle sub-plot is angular rate vs time
% The botom sub-plot is angular acceleration vs time
% All axis should be clearly labeled (including units)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Make plots

% Set up the figure:
figure(1030); clf;
subplot(211); hold on;grid on;grid minor;
plot(tOde45,zOde45(1,:),'r-','LineWidth',2)
plot(tGrid,zEuler(1,:),'b-','LineWidth',2)
xlabel('time(s)');
ylabel('angle(rad)');
title('Pendulum Simulation');
legend('ode45', 'euler');

subplot(212);hold on;grid on;grid minor;
plot(tOde45,zOde45(2,:),'r-','LineWidth',2)
plot(tGrid, zEuler(2,:),'b-','LineWidth',2)
xlabel('time (s)');
ylabel('rate (rad/s)');
legend('ode45', 'euler');
end
