function [x_opt,f_opt,x_kept,f_kept] = DSA(MyFunc,lb,ub,DSA_options)

% -----------------------------------------------------------------
% Function implementing the differential simulated annealing (DSA)
% algorithm for nonlinear least squares fitting
% -----------------------------------------------------------------
% Input parameters
% -----------------------------------------------------------------
% 1. MyFunc - function handle for the residues to be minimized. It takes the
% form res = MyFunc(x), in which x is the vector storing variables to be
% optimized, res is a vector storing the residues. The objective function
% value is the sum of squares of elements in res.
% 2. lb, ub - vectors storing lower and upper bounds for elements of x
% 3. DSA_options - parameters used in the optimization process. It is a
% structure with fields start_beta, end_beta, cool_rate, local_step,
% mc_length, expected_fobj, and visual_output. If this argument is not provided,
% default parameters will be used.
% -----------------------------------------------------------------
% Outputs
% -----------------------------------------------------------------
% 1. x_opt - optimal solution found
% 2. f_opt - optimal objective function value found
% 3. x_kept - solutions kept for satisfying fitting performance (defined as
% objective function value < 10*EXPECTED_FOBJ
% 4. f_kept - objective function values corresponding to x_kept
% -----------------------------------------------------------------

% Define global parameters
global f_now fopt_now xopt_now counts_accepted counts_all nr

% Set up parameters for optimization
if nargin == 3
    % No input for DSA_options, use the default parameters
    START_BETA=1e-10; 					% Initial temperature
    END_BETA=1e10;    					% Maximal temperature
    COOLING_RATE=1.03;					% Cooling rate
    LOCAL_STEP=0.05;					% Step size for local perturbation
    MC_LENGTH=fix(50*(length(lb)^0.5)); % Markov chain length
    EXPECTED_FOBJ=1e-6;					% Expected value of objective function
    VISUAL_OUTPUT=1; 					% Show convergence curve (1) or textual output only (0)
elseif nargin == 4
    % Use parameters specified in DSA_options
    if isfield(DSA_options,'start_beta')
        START_BETA=DSA_options.start_beta;
    else
        START_BETA=1e-10;
    end
    if isfield(DSA_options,'end_beta')
        END_BETA=DSA_options.end_beta;
    else
        END_BETA=1e10;
    end
    if isfield(DSA_options,'cooling_rate')
        COOLING_RATE=DSA_options.cooling_rate;
    else
        COOLING_RATE=1.03;
    end
    if isfield(DSA_options,'local_step')
        LOCAL_STEP=DSA_options.local_step;
    else
        LOCAL_STEP=0.05;
    end
    if isfield(DSA_options,'mc_length')
        MC_LENGTH=DSA_options.mc_length;
    else
        MC_LENGTH=500;
    end
    if isfield(DSA_options,'expected_fobj')
        EXPECTED_FOBJ=DSA_options.expected_fobj;
    else
        EXPECTED_FOBJ=0;
    end
    if isfield(DSA_options,'visual_output')
        VISUAL_OUTPUT=DSA_options.visual_output;
    else
        VISUAL_OUTPUT=1;
    end
else
    fprintf('Incorrect number of input parameters\n');
    x_opt = NaN;
    f_opt = NaN;
    return;
end

% Generate output
fprintf('----------------------------\n');
fprintf('Global optimization with DSA\n');
fprintf('----------------------------\n');
fprintf('Initial beta=%.2e\n',START_BETA);
fprintf('Maximal beta=%.2e\n',END_BETA);
fprintf('Cooling rate=%f\n',COOLING_RATE);
fprintf('Inital step size=%f\n',LOCAL_STEP);
fprintf('Markov chain length=%d\n',MC_LENGTH);
fprintf('Expected optimal objective function=%.2e\n',EXPECTED_FOBJ);
fprintf('----------------------------\n');

% Set up global variables
np=length(lb);
x = (lb+ub)/2; %Start from the center of feasible region
f_now = sumsqr(MyFunc(x));
fopt_now = f_now;
xopt_now = x;
counts_accepted = zeros(np,1);
counts_all = zeros(np,1);
nr = length(MyFunc(x));
local_steps = LOCAL_STEP*ones(np,1);
beta = START_BETA;


% Visualize output
if VISUAL_OUTPUT
    h = animatedline;
    axis([START_BETA END_BETA EXPECTED_FOBJ 1.2*fopt_now]);
    title('Convergence of DSA');
    xlabel('Temperature \beta');ylabel('Best objective function value found');
    set(gca,'XScale','log','YScale','log');
end

% Iterations for the main phase of optimization
n_func_eval=0;
i_x_kept=0;
x_kept=zeros(1000,np); %Pre-allocate space for x_kept and f_kept
f_kept=zeros(1000,1);
while beta < END_BETA && fopt_now > EXPECTED_FOBJ
    for n=1:MC_LENGTH
        n_func_eval = n_func_eval + 1;
        if mod(n_func_eval,500) == 1            
            [eig_vecs,eig_vals] = getDirections(MyFunc,x);
            if n_func_eval > 1
                local_steps = adjustStep(local_steps,counts_accepted,counts_all);
                counts_accepted = zeros(np,1);
                counts_all = zeros(np,1);
            end
        end
        x = MonteCarloStep(MyFunc,x,lb,ub,local_steps,eig_vecs,eig_vals,beta);
        if f_now < 100*EXPECTED_FOBJ
            i_x_kept = i_x_kept+1;
            x_kept(i_x_kept,:) = x';
            f_kept(i_x_kept) = f_now;
        end
    end
    if VISUAL_OUTPUT
    addpoints(h,beta,fopt_now);
    drawnow;
    else
        fprintf('Beta = %.2e\n',beta);
        fprintf('Optimal objective function value = %.2e\n',fopt_now);
    end
    beta = beta*COOLING_RATE;
end
x_opt = xopt_now;
f_opt = fopt_now;
x_kept(f_kept==0,:)=[];
f_kept(f_kept==0)=[];
[x_kept,ia,~]=unique(x_kept,'row');
f_kept=f_kept(ia);
end

function H=approxHamiltonian(MyFunc,x)
% Calculate approximate Hamiltonian based on Jacobian of the residues with
% respect to the variables to be optimized
global nr
h=0.01;
np=length(x);
J=zeros(nr,np);
for i=1:np
    xr=x;xr(i)=xr(i)+h;xl=x;xl(i)=xl(i)-h;
    J(:,i)=(MyFunc(xr)-MyFunc(xl))/(2*h);
end
H=2*J'*J;
end

function [eig_vecs,eig_vals]=getDirections(MyFunc,x)
% Calculate eigenvalues and eigenvectors of the approximate Hamiltonian
H=approxHamiltonian(MyFunc,x);
[eig_vecs,D]=eig(H);
eig_vals=diag(D);
end

function r=getMaxStep(x0,lb,ub,dx)
% Compute the range of step size for the perturbation x1=x0+k*dx, s.t. lb<=x1<=ub
d_bounds=[lb-x0 ub-x0]./dx;
d_lb=min(d_bounds')';
d_ub=max(d_bounds')';
r=[max(d_lb) min(d_ub)];
end

function x1=Perturbation(x0,lb,ub,step_local,dx,label_perturb)
%Perturb the current solution to generate a new candidate solution
% x0 - current solution
% lb,ub - lower and upper bounds for solutions
% step_local - current step size for local perturbation
% dx - direction for perturbation
% label_perturb - type of perturbation: 0 - global perturbation; 1 - local perturbation
rglobal=getMaxStep(x0,lb,ub,dx);
if label_perturb == 1
    rlocal=[max(-step_local,rglobal(1)) min(step_local,rglobal(2))];
    k=rand()*(rlocal(2)-rlocal(1))+rlocal(1);
    x1=x0+k*dx;
else
    k=rand()*(rglobal(2)-rglobal(1))+rglobal(1);
    x1=x0+k*dx;
end
end

function r = calcRatio(x)
%Calculate the coefficients for step size adjustment
c=3;
if x<0.4
    r=1/(1+c*(0.4-x)/0.4);
elseif x>0.6
    r=1+c*(x-0.6)/0.4;
else
    r=1;
end
end

function s = adjustStep(s0,cnts_accepted,cnts_all)
% Adjust step size for local perturbation based on the fraction of accepted perturbations on this direction
frac_accepted = cnts_accepted./cnts_all;
rvec = arrayfun(@calcRatio,frac_accepted);
rvec(cnts_all == 0)=1;
s = s0.*rvec;
end

function x1=MonteCarloStep(MyFunc,x0,lb,ub,local_steps,eig_vecs,eig_vals,beta)
% Execute one step of Metropolis sampling from solution x0
global f_now fopt_now xopt_now counts_accepted counts_all
np=length(x0);
idx_perturb = randi(np); %Pick the idx_perturb-th eigenvector as the direction for perturbation
if eig_vals(idx_perturb)<1000/beta
    label_perturb = 0;
else
    label_perturb = 1;
    counts_all(idx_perturb) = counts_all(idx_perturb)+1;
end
x1=Perturbation(x0,lb,ub,local_steps(idx_perturb),eig_vecs(:,idx_perturb),label_perturb);
f0=f_now;
%MyFunc(x1)
if isreal(sum(x1))
    new_residues = MyFunc(x1);
    f1 = sumsqr(new_residues);
    p_accept=exp((f0-f1)*beta);
else
    f1 = 1e20;
    p_accept = 0;
end
%[f0 f1]

if rand() < p_accept
    f_now=f1;
    counts_accepted(idx_perturb) = counts_accepted(idx_perturb)+1;
    if f1 < fopt_now
        fopt_now = f1;
        xopt_now = x1;
    end
else
    x1=x0;
end
end