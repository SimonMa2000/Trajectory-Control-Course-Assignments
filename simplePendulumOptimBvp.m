function soln = simplePendulumOptimBvp(config, param, nlpOpt)
% Probably 10 hours.
% soln = simplePendulumOptimBvp(config, param, nlpOpt)
%
% Compute the minimum-torque solution to move the simple pendulum between
% two prescribed states in a specified duration. The transcription of the
% trajectory optimization problem is performed using Euler's method
% and multiple shooting with one step per segment. The duration of each
% step is the same (uniform time grid).
%
% INPUTS:
%   config = struct = configuration options for the trajectory optimization
%       .nStep = scalar = number of simulation steps
%       .beginState = [2, 1] = initial [angle; rate]
%       .finalState = [2, 1] = final [angle; rate]
%       .duration = scalar = duration of the trajectory
%   param = struct = parameters of the simple pendulum
%     .freq = scalar = undamped natural frequency squared
%                    = (gravity / length) for a point mass pendulum
%     .damp = scalar = normalized linear viscous friction term
%   nlpOpt = solver options object, created by:
%        >> nlpOpt = optimset('fmincon')
%       Useful options (set using "." operator, eg: nlpOpt.Display = 'off')
%           --> Display = {'iter', 'final', 'off'}
%           --> OptimalityTolerance = {1e-3, 1e-6}
%           --> ConstraintTolerance = {1e-3, 1e-5, 1e-10}
%
% OUTPUTS:
%   soln = struct = solution data
%    .grid = struct = values of the trajectory at the grid points
%       .time = [1, nStep+1] = knot points
%       .state = [2, nStep+1] = state at the knot points
%       .control = [1, nStep] = control over each step (constant)
%   .info = information about the optimization run
%       .nlpTime = time (seconds) spent in fmincon
%       .exitFlag = fmincon exit flag
%       .objVal = value of the objective function
%       .[all fields in the fmincon "output" struct]
%
% NOTES:
%
%   Minimize: integral of actuator torque along the entire trajectory.
%   J = \int_0^T u^2(t) dt
%
%   Subject to: system dynamics and boundary conditions
%

%%%%%
nSegment = config.nStep;
nSubstep = 1;
nControl = 1;
nState = 2;
if (nSegment<2) 
    error('Need at least 2 segments');
end
T = config.duration;
% set up grid
tGrid = linspace(0,T,nSegment+1); 
nGrid = length(tGrid);
%%% initial guesses
% Size of  simple initial guesses
% guess.z = [1*nGrid]
% guess.u = [1*nGrid-1]
% Tried different initial guess - may get better result

% Heuristic guess 1:
% guess.z = [linspace(0,pi,tGrid);
% Heuristic guess 2: 
% guess.z = [ -pi/config.duration.*tGrid.*sin(param.freq*tGrid); 
%             -pi/config.duration.*tGrid.*sin(param.freq*tGrid-pi/2);];
% Heuristic guess 3: 
guess.z = [ -pi*sin(pi/2*config.duration.*tGrid).*sin(param.freq*tGrid); 
          -pi*sin(pi/2*config.duration.*tGrid).*sin(param.freq*tGrid-pi/2);];
% keep guess.u as zero vector (zero torque).
% The control on the last time step would not
% be decision variable in multiple shooting method
guess.u = zeros(1,nGrid);
guess.uDec = guess.u(1:end-1);
%%% nonlinear programming solver
problem.x0 = pack(guess.z, guess.uDec, nGrid, nState, nControl, nSubstep);
problem.lb = []; problem.ub = [];
problem.Aineq = []; problem.bineq = [];
problem.Aeq = []; problem.beq = [];
problem.objective = @(decVar)objective(decVar,nGrid,tGrid);
problem.nonlcon = @(decVar)nonLinCst(decVar, tGrid(1:end-1), config, param); %    
problem.solver = 'fmincon';
problem.options = nlpOpt; 
[zuSoln, fVal, ~] = fmincon(problem);
%%% Store the trajectory in a nice format: 
soln.grid.time = tGrid;
% soln.grid.state(1,:) = zSoln(1,:);
% soln.grid.state(2,:) = zSoln(2,:);
% soln.grid.control = uSoln, last point control deleted!
[soln.grid.state,soln.grid.control] = unpack(zuSoln, nGrid, nState, nControl, nSubstep);
soln.info.objVal = fVal; 
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
function cost = objective(decVar,nGrid,tGrid)      
             [~,u] = unpack(decVar, nGrid);
             h = tGrid(2)-tGrid(1);
             cost = h * sum(u.^2);
end
function [C, Ceq] = nonLinCst(decVar, tGrid, config, param); % 
    nState = 2; nControl = 1; 
    nGrid = config.nStep+1 ; nSubstep = 1;
    [z, uDec] = unpack(decVar, nGrid, nState, nControl, nSubstep); % delete time as decVar z
    % let the last term of uDec hold for one more step
    u = [uDec,uDec(end)];
    % hard code euler method 
    h = tGrid(2)-tGrid(1);
    zNext = z + h * simplePendulumDynamics(z, u ,param);
    % Index Info:  z(nState,nGrid,nSubStep)
    % z is result of Euler method
    %%% Boundary Value Constraints:
    InitCst = z(:,1) - config.beginState; %Initial Position 
    FinalCst = z(:,end) - config.finalState; %Final Position
    %%% Defect Constraints:
    zEnd = zNext(:, 1:(end-1));  % States at the end of a segment
    zStart = z(:, 2:end);   % States at the beginning of a segment
    Defect = reshape(zStart-zEnd, nState*(nGrid-1),1);
    C = [];  % No inequality constraints
    Ceq = [InitCst;Defect;FinalCst]; % Boundary Condition
end



% z [nState, nGrid*nSubstep] state vector
% u [nControl,(nGrid -1)*nSubstep] control vector
% decVar [ nState*nGrid*nSubstep; nControl*(nGrid*nSubstep-1)]
% 
function [z,u]= unpack(decVar, nGrid, nState, nControl, nSubstep)
% z0 decVar in human readable matrix form
% z0 = [[q;w];u](mat) where q,w,u are row vector
% decVar is [q,w,u] row vector
if (nargin == 2)
    nState = 2;
    nControl = 1;
    nSubstep = 1;
end
% add the last term of u for batch process
decVar = [decVar; decVar(end)];
% batch process
 z0_trans = reshape(decVar,(nGrid)*nSubstep,nState+nControl);
 z0 = z0_trans';
 z = z0(1:nState,:);
% (end-1) to delete the last control term
u = z0(nState+1:end,1:end-1);
end

function decVar = pack(z, u, nGrid, nState, nControl, nSubstep)
% decVar column vector [(nControl+nState)*nGrid*nSubstep,1]
if (nargin == 3)
    nState = 2;
    nControl = 1;
    nSubstep = 1;
end   
% add last control term
u = [u, u(end)];
decVar = reshape([z;u]', (nControl+nState)*nGrid*nSubstep,1);
% delete last control term
decVar(end)=[];
end