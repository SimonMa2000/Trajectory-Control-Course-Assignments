function dState = LorenzDynamics(~, state)
%  dState = LorenzDynamics(time, state)
%
% This function computes the dynamics for the famous 'Lorenz Attractor'
% See the matlab function lorenz.m for a neat demo.
%
% INPUTS:
%   time = [1, nTime] = row vector of query times (unused)
%   state = [3, nTime] = current state of the system
%
% OUTPUTS:
%   dState = [3, nTime] = time-derivative at current state
%
% NOTES:
%   https://en.wikipedia.org/wiki/Lorenz_system
%
% dx = sigma * (y - x)
% dy = x * (rho - z) - y
% dz = x * y - beta * z
%

%% Implement this function
% unpack State
x = state(1, :);
y = state(2, :);
z = state(3, :);

sigma = 10;
rho = 29;
beta = 8/3;

% compute derivatives
dx = sigma * (y - x);
dy = x * (rho - z) - y;
dz = x * y - beta * z;

% pack dState
dState = [dx; dy; dz];

end
