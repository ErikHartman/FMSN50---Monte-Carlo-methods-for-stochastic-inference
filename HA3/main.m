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
rho=0.05;
burn_in = 10;
M=1000;
N=burn_in+M; % Nbr iterations
thetas = zeros(N+1,1); % Initialize vectors of theta
lambdas = zeros(N+1,d); % initialize matrix of lambda
breakpoints=create_breakpoints(d,T);
n_tau= nbr_event_between_breakpoints(T, breakpoints);
thetas(1)=init_theta(rho); % Init theta
lambdas(1,:)=init_lambdas(thetas(1), d); % Init lambda
tries = zeros(N,1);
matrix_breakpoints = zeros(N, d+1);
accepted_vector = zeros(N,1);
% Main algorithm, k steps
for k=1:N % Main loop
    % Start with MH! Advice from Isabella
    [accepted, new_breakpoints] = MH_algorithm(lambdas(k,:), breakpoints, T,rho);
    accepted_vector(k) = accepted;
    matrix_breakpoints(k,:) = new_breakpoints;
    breakpoints = new_breakpoints;
    % breakpoints(:,i+1) = new_breakpoints;
    % Draw theta Gibbs
    thetas(k+1)=sample_theta(d, rho, lambdas(k,:));
    % Draw lambdas Gibbs
    lambdas(k+1,:) = sample_lambdas(T, thetas(k+1), breakpoints, d);
end
theta_est=thetas(burn_in:end);
lambdas_est = lambdas(burn_in:end, :);
breakpoints_est=matrix_breakpoints(burn_in:end,:);
%% Task 2c
% Plot and save images for d=2,3,4,5, burn_in=1 for lambda
d=5;
rho=0.05;
burn_in=1;
M=10000;
[theta_est, lambdas_est, breakpoints_est]=main_script(d,T,rho,burn_in, M, rho);
xlabel("Iterations")
ylabel("Lambda")
hold on
for i=1:d
    str=['Lambda: ' num2str(i)]
    plot(lambdas_est(:, i),'DisplayName',str)
end
title(['Nbr breakpoints: ' num2str(d-1)])
hold off
legend show
%% Plots
d=5;
rho=0.05;
burn_in=1;
M=10000;
[theta_est, lambdas_est, breakpoints_est]=main_script(d,T,rho,burn_in, M, rho);
figure;
hold on
histogram(T, 50)
ylabel("Frequency")
for i=1:d+1
    %xline(breakpoints_est(1,i), 'r--');

    xline(breakpoints_est(end,i));
end
hold off

%% Plot breakpoint trajectory, be careful abou notation
d=5;
rho=0.05;
burn_in=1;
M=10000;
[theta_est, lambdas_est, breakpoints_est, ~]=main_script(d,T,rho,burn_in, M, rho);
figure;
hold on
for k = 2:d
    str = ['Breakpoint: ' num2str(k-1)];
    plot(breakpoints_est(:,k), 'DisplayName', str);
end
title(['Number of breakpoints: ' num2str(d-1)])
hold off
legend show
%% Plot histogram of parameters
figure;
hold on;
for k = 2:d
    histogram(lambdas(:,k));
    histogram(thetas);    
end

hold off

%% 2d. Try different psi

d=5;
rho=0.05;
burn_in=100;
M=1000;
len = 50;
psi_x=1:1:50;
mean_theta = zeros(len,1);
var_theta = zeros(len,1)
mean_lambda = zeros(len,d);
var_lambda = zeros(len,d);
for rho=1:len
    [theta_est, lambda_est,~, ~ ]=main_script(d,T,rho,burn_in, M, rho);
    mean_theta(rho)=mean(theta_est);
    var_theta(rho) = var(theta_est);
    mean_lambda(rho, :) = mean(lambda_est)';
    var_lambda(rho, :) = var(lambda_est);
end
%%
figure;
plot(psi_x, mean_theta, 'bo')
xlabel("\psi")	
ylabel("Mean \theta")	
set(gca,'FontSize',14)

figure;
plot(psi_x, var_theta, 'bo')
xlabel("\psi")	
ylabel("Variance \theta")	
set(gca,'FontSize',14)
%%
figure;
hold on
for i=1:d
    str=['\lambda' num2str(i)];
    plot(psi_x,mean_lambda(:,i), '.', 'MarkerSize', 14, 'DisplayName', str)
end
hold off
xlabel("\psi")	
ylabel("Mean \lambda")	
set(gca,'FontSize',14)
legend show

%%
figure;
hold on
for i=1:d
    str=['\lambda' num2str(i)];
    plot(psi_x,var_lambda(:,i), '.', 'MarkerSize', 14, 'DisplayName', str)
end
hold off
xlabel("\psi")	
ylabel("Variance \lambda")	
set(gca,'FontSize',14)
legend show

%% 2e change rho and plot acceptance probability
d=5;
burn_in=100;
M=1000;
rho_x=0.001:0.001:0.1;
len = length(rho_x);
mean_acceptance = zeros(len,1);
var_acceptance = zeros(len,1);

for rho=1:len
    [theta_est, lambda_est, ~, accepted_vector_est]=main_script(d,T,rho,burn_in, M, rho);
    mean_acceptance(rho)=mean(accepted_vector_est);
    var_acceptance(rho) = var(accepted_vector_est);

end

%%
figure;
plot(rho_x, mean_acceptance, 'bo')
xlabel("\rho")
ylabel("Mean ACR")
set(gca,'FontSize',14)

figure;
plot(rho_x, var_acceptance, 'bo')
xlabel("\rho")
ylabel("Variance ACR")
set(gca,'FontSize',14)


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
    init_theta = gamrnd(2,1/psi);
end

% Initialize lambdas by prior distribution, returns a full vector.
function init_lambdas=init_lambdas(theta, nbr_breakpoints)
    init_lambdas = gamrnd(2,1/theta, nbr_breakpoints, 1);

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
    n_tau = nbr_event_between_breakpoints(tau, breakpoints);
    for i=1:nbr_breakpoints
        time_difference(i) = breakpoints(i+1)-breakpoints(i);
        lambdas(i) = gamrnd(2+n_tau(i), 1/(theta+time_difference(i)))'; % time_difference= t_i+1-t_i
    end
end

function proposed_breakpoint=random_walk_proposal(breakpoints, rho, index)
    % Contruct R
    R=rho*(breakpoints(index+1)-breakpoints(index-1));
    epsilon = unifrnd(-R,R);
    proposed_breakpoint = breakpoints(index)+epsilon;

end

% Function calculate transition probability and compare with 1. Here we use
% symmetric proposal kernel, meaning we sample from 1^f(z)/f(x) (z=t*, x=t)
% This is in our case 1^(f(t*|tau, lambda, theta)/f(t|tau, lambda, theta))
% From exercise a this simplifies to:
% 1^(f(t*)*f(tau|lambda, t*)/f(t)*f(tau|lambda, t)) <- These we can
% sample!!
function acceptance_probability = calculate_accaptance_probability(lambdas_k,proposed_breakpoints, breakpoints, tau)
    % For calculating f(t*) f(tau|lambda, t*), we need: 
    % Number of disasters for proposed, lambdas for iteration k and
    % proposed breakpoints:
    if issorted(proposed_breakpoints)
        n_proposed = nbr_event_between_breakpoints(tau, proposed_breakpoints);
        n_k = nbr_event_between_breakpoints(tau, breakpoints);
        diff_proposed = proposed_breakpoints(2:end)-proposed_breakpoints(1:end-1);  % t_i+1 - t_i
        diff_breakpoints = breakpoints(2:end)-breakpoints(1:end-1);

        % Log since otherwise big and messy numbers
        % logaritmera analytiskt
        log_r_t_star = sum(log(diff_proposed))-sum(lambdas_k'.*diff_proposed)+sum(log(lambdas_k').*n_proposed);
        log_r_t = sum(log(diff_breakpoints))-sum(lambdas_k'.*diff_breakpoints)+sum(log(lambdas_k').*n_k);
        acceptance_probability = min(1,exp(log_r_t_star-log_r_t));
    else
        acceptance_probability = 0;
    end
end

% Function return new sampled breakpoints for each iteration, see slide 5
% Lecture 10. breakpoints is the breakpoints input for iteration k.
function [accepted, ret] = MH_algorithm(lambdas_k, breakpoints, tau, rho)
    % Start by copying the old breakpoints to next step since only a
    % fraction of them will update.
    % Main for loop, go through all breakpoints and propose new ones.
    for i=2:length(breakpoints)-1
      suggested_breakpoints = breakpoints;
      % Propose a random walk for each breakpoint
      proposed_breakpoint=random_walk_proposal(suggested_breakpoints, rho, i);
      suggested_breakpoints(i) = proposed_breakpoint;
      % Calculate probability of accepting the proposed breakpoint and add
      % to array.
      acceptance_probability = calculate_accaptance_probability(lambdas_k,suggested_breakpoints, breakpoints, tau);
      % Calculate acceptance by draw from uniform and compare with
      % acceptance rate.
       % If accepted, update breakpoint to the proposed one
        accepted = 0;
      alpha = unifrnd(0,1);
      if alpha <= acceptance_probability
          breakpoints = suggested_breakpoints;
          accepted = 1;
      end
    end
    ret = breakpoints;
end

% Main script, d = mbr breakepoints, T=tau, psi=hyperparameter,
% burn_in=iterations to rm from estimate, M= actrual estimate iterations.
function [theta_est, lambdas_est, breakpoints_est, accepted_vector_est]=main_script(d,T,psi,burn_in, M, rho)
    N=burn_in+M; % Nbr iterations
    thetas = zeros(N+1,1); % Initialize vectors of theta
    lambdas = zeros(N+1,d); % initialize matrix of lambda
    breakpoints=create_breakpoints(d,T);
    n_tau= nbr_event_between_breakpoints(T, breakpoints);
    thetas(1)=init_theta(psi); % Init theta
    lambdas(1,:)=init_lambdas(thetas(1), d); % Init lambda
    matrix_breakpoints = zeros(N, d+1);
    accepted_vector = zeros(N,1);
    % Main algorithm, k steps
    for k=1:N % Main loop
        % Start with MH! Advice from Isabella
        [accepted, new_breakpoints] = MH_algorithm(lambdas(k,:), breakpoints, T, rho);
        accepted_vector(k) = accepted;
        matrix_breakpoints(k,:) = new_breakpoints;
        breakpoints = new_breakpoints;
        % breakpoints(:,i+1) = new_breakpoints;
        % Draw theta Gibbs
        thetas(k+1)=sample_theta(d, psi, lambdas(k,:));
        % Draw lambdas Gibbs
        lambdas(k+1,:) = sample_lambdas(T, thetas(k+1), breakpoints, d);
    end
    accepted_vector_est = accepted_vector(burn_in:end);
    theta_est=thetas(burn_in:end);
    lambdas_est = lambdas(burn_in:end, :);
    breakpoints_est=matrix_breakpoints(burn_in:end,:);
end
