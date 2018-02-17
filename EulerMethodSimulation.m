function [tGrid, zGrid] = EulerMethodSimulation(dynFun, tSpan, zInit, hMax)
% [tGrid, zGrid] = EulerMethodSimulation(dynFun, tSpan, zInit, hMax)
%
% Simulate a dynamical system over a set time span using Euler's method,
% given the initial state and the maximum time step. Use a uniform time grid.
%
% INPUTS:
%    dynFun = a function handle:  dz = dynFun(t, z)
%        IN:  t = [1, nTime] = row vector of time
%        IN:  z = [nState, nTime] = matrix of states corresponding to each time
%        OUT: dz = [nState, nTime] = time-derivative of the state at each point
%    tSpan = [startTime, finalTime] = [1, 2] = time span
%    zInit = [nState, 1] = state of the system at start time
%    hMax = scalar = maximum time step to use for the uniform time grid.
%
% OUTPUTS:
%    tGrid = [1, nTime] = time grid on which the integration was performed
%    zGrid = [nState, nTime] = state at each point in tGrid
%
%
% Implement this function
% the following code are modified from DEMO in problem 1,Lecture 3
nTime=size(tSpan,2);
nState=size(zInit,1);
zGrid=zeros(nState,nTime);
zGrid(:,1)=zInit;
% shift tGrid to initial tGrid=0; 
tGrid=tSpan-tSpan(1)*ones(1,nTime);
for iGrid=1:(nTime-1)
    % local variable used in Euler iteration,nTimes-1 times iterstion
    tPrev=tSpan(:,iGrid);
    zPrev=zGrid(:,iGrid);
    
    % compute dynamic function
    dzPrev=dynFun(tPrev,zPrev);
    
    % Euler's iteration
    zGrid(:, iGrid+1) = zPrev + hMax * dzPrev; 
end
end
