%% Part 1
load('coal_mine_disasters.mat')

% Plot the data
%histogram (T)

% Define initial gamma(2,psi) to sample our theta, start by setting psi=1.
% f(theta) in the task
n=2000
gam_theta = @(psi,n) gamrnd(2,psi,n,1)
thetas = gam_theta(2,n);

% f(lambda|theta)
gam_lambda_theta = @(theta) gamrnd(2,theta,n,1)
lambdas_given_theta = gam_lambda_theta(thetas);

% f(lambda|theta)*f(theta)
product=lambdas_given_theta.*thetas;

histogram(product) % Gamma, we need to find the paramaters for each of these!


%% Main script
load('coal_mine_disasters.mat')

d=2; % Nbr breakpoints
psi=2;
N=1; % Nbr iterations
thetas = zeros(N+1,1); % Initialize vectors of theta
lambdas = zeros(N+1,d); % initialize matrix of lambda
breakpoints=create_breakpoints(d,T);
n_tau= nbr_event_between_breakpoints(T, breakpoints);
thetas(1)=init_theta(psi); % Init theta
lambdas(1,:)=init_lambdas(thetas(1), d); % Init lambda

% Main algorithm
for i=1:N % Main loop
    % Draw theta Gibbs
    thetas(i+1)=sample_theta(d, psi, lambdas(i,:));
    % Draw lambdas Gibbs
    lambdas(i+1,:) = sample_lambdas(T, thetas(i+1), breakpoints, d);
    
    % Update breakpoints, MH
        % Loop through all breakpoints.
        % Generate random walk proposal
        % Calculate acception probability
        % Reject or accept the new breakpoint
        % Update breakpoint.
    
end


%% Useful functions
% Calculates nbr of event between breakpoints returns the full vector for
% each breakpoint.


function t_out=create_breakpoints(d,tau)
    t_out = linspace(tau(1), tau(end), d+1)';
end

function n_tau = nbr_event_between_breakpoints(tau, breakpoints)
    for i=1:length(breakpoints)-1
        n_tau(i) = sum(tau <= breakpoints(i+1) & tau>=breakpoints(i)); % Larger than or equal larger than?
    end
end

% Initialize theta by prior distribution:
function init_theta = init_theta(psi)
    init_theta = gamrnd(2,psi);
end

% Initialize lambdas by prior distribution, returns a full vector.
function init_lambdas=init_lambdas(theta, nbr_breakpoints)
    init_lambdas = gamrnd(2,theta, nbr_breakpoints, 1);

end

% Function to sample from conditional theta based on posterior and prior
% In other words: f(theta|lambda,t,tau)= f(theta)*f(lambda|theta)
function theta = sample_theta(nbr_breakpoints, psi, lambda)
    % The conditional theta is a gamma distributed, see report. Thus we
    % sample:
    theta = gamrnd(2*nbr_breakpoints + 2, psi + sum(lambda));
end

% Function to sample from conditional lambda based on posterior
% distributions. In other words:
% f(lambda|t,theta,tau)=f(tau|lambda,t)*f(lambda|theta)

function lambdas = sample_lambdas(tau, theta, breakpoints, nbr_breakpoints)
    n_tau = nbr_event_between_breakpoints(tau, breakpoints)
    for i=1:nbr_breakpoints
        time_difference(i) = breakpoints(i+1)-breakpoints(i);
        lambdas(i) = gamrnd(2+n_tau(i), theta+time_difference(i))'; % time_difference= t_i+1-t_i
    end
end

