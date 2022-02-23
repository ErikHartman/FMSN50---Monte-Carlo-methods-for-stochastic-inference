%% Part 1
load('coal_mine_disasters.mat')

% Plot the data
%histogram (T)

% Define initial gamma(2,psi) to sample our theta, start by setting psi=1.
% f(theta) in the task
n=200
thetas = gam_theta(2,n);
gam_theta = @(psi,n) gamrnd(2,psi,n,1)
% f(lambda|theta)
gam_lambda_theta = @(theta) gamrnd(2,theta,n,1)
lambdas_given_theta = gam_lambda_theta(thetas);

% f(lambda|theta)*f(theta)
product=lambdas_given_theta.*thetas;

histogram(product) % Exponential?

%%  f(tau|lambda,t)*f(lambda|theta)

% f(lambda|theta) = product
% f(tau|lambda,t)

n_i = @(tau,i) sum(tau(tau>i & tau<i+1));


tau_given_lambda_theta = @(t,lambda, tau) exp(-sum(lambda.*(t(2:end)-t(1:end-1))))%.*cumprod(lambda)
n=3 % Three breakpoits
t_out=create_breakpoints(n,T)
thetas = gam_theta(2,n);
gam_lambda_theta = @(theta) gamrnd(2,theta,n,1)
lambdas = gam_lambda_theta(thetas);

tau_given_lambda_theta(t_out, lambdas, T)
nbr_event_between_breakpoints(1, T, t_out)
% Function to create breakpoints, e.g evenly spaces points in time.
function t_out=create_breakpoints(d,tau)
    t_out = linspace(tau(1),tau(end),d)
end

% Calculates nbr of event between breakpoint t, i and i+1
function n_i = nbr_event_between_breakpoints(i, tau, t)
    start_time = t(i)
    end_time = t(i+1)
    n_i = sum(tau > start_time & tau < end_time)
end


