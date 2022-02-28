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
N=3; % Nbr iterations
thetas = zeros(N+1,1); % Initialize vectors of theta
lambdas = zeros(N+1,d); % initialize matrix of lambda
breakpoints=create_breakpoints(d,T);
n_tau= nbr_event_between_breakpoints(T, breakpoints);
thetas(1)=init_theta(psi); % Init theta
lambdas(1,:)=init_lambdas(thetas(1), d); % Init lambda
tries = zeros(N,1);

% Main algorithm
for i=1:N % Main loop
    % Start with MH!

   % breakpoints(:,i+1) = new_breakpoints;
    % Draw theta Gibbs
    thetas(i+1)=sample_theta(d, psi, lambdas(i,:));
    % Draw lambdas Gibbs
    lambdas(i+1,:) = sample_lambdas(T, thetas(i+1), breakpoints, d);
    % Update breakpoints, MH
    nr_current_tries = 0;
    accepted = 0;
    while ~accepted % oops! This never ends! Something is wrong with lambda or something
     [new_breakpoints, accepted] = MH(lambdas(i,:), breakpoints(:,i), T);
      tries(i)  = tries(i) + 1; % want around 30% acceptance rate
    end

end


%% Plot histogram with breakpoints
figure;
hold on
histogram(T)
title('Disasters')
xlabel('t')
for i = 1:length(breakpoints)
    xline(breakpoints(i));
end
hold off

%% Plot breakpoint trajectory
figure;
hold on
for i = 1:d
    plot(breakpoints(i,:));
end
%% Plot histogram of parameters
figure;
hold on;
histogram(lambdas);
histogram(thetas);
hold off

%% Useful functions
% Calculates nbr of event between breakpoints returns the full vector for
% each breakpoint.


function t_out=create_breakpoints(d,tau)
    t_out = linspace(tau(1), tau(end), d+1)';
end

function n_tau = nbr_event_between_breakpoints(tau, breakpoints)
    n_tau = zeros(length(breakpoints)-1,1);
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
    theta = gamrnd(2*nbr_breakpoints + 2, 1/(psi + sum(lambda)));
end

% Function to sample from conditional lambda based on posterior
% distributions. In other words:
% f(lambda|t,theta,tau)=f(tau|lambda,t)*f(lambda|theta)

function lambdas = sample_lambdas(tau, theta, breakpoints, nbr_breakpoints)
    n_tau = nbr_event_between_breakpoints(tau, breakpoints)
    for i=1:nbr_breakpoints
        time_difference(i) = breakpoints(i+1)-breakpoints(i);
        lambdas(i) = gamrnd(2+n_tau(i), 1/(theta+time_difference(i)))'; % time_difference= t_i+1-t_i
    end
end

% Function to sample from conditional t for MH-agorithm
% f(t|theta, lambda, T) 
function fX = f(lambda,t,T)
    startpoints = t(1:end-1);
    endpoints = t(2:end);
    accidents = zeros(length(t)-1, 1);
    for j = 1:length(startpoints)
        accidents(j) = sum(length(T(T > startpoints(j) & T < endpoints(j))));
    end
    fX = exp(sum(log(lambda).*accidents + log(t(2:end)-t(1:end-1)) - lambda.*(t(2:end)-t(1:end-1))));
end

% Lecture 10 slide 5 har pseudokod
% Tog från Bolin
function [breakpoints, accepted] = MH(lambda, breakpoints, T)
    accepted = 0;
    % Selecting one t_i at random
    i = randi(length(breakpoints) - 2) + 1;
    % Calculate the maximum step size based on the interval width
    rho = 0.055;
    R = rho*(breakpoints(i+1)-breakpoints(i-1)); % random walk proposal all at once
    X = -inf;
    % Keep drawing until we get an X that is inside the interval
    while(X < breakpoints(i-1) || X > breakpoints(i+1))
        X = breakpoints(i) + R*(2*rand(1)-1);
    end
    new_br = breakpoints;
    new_br(i) = X;
    % Accept?
    if rand(1) <= f(lambda, new_br, T)/f(lambda, breakpoints, T)
        breakpoints = new_br;
        accepted = 1;
    end
end

