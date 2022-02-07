

%% Parameters
lambda_vec = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; % A
k_vec = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0]; % B
load powercurve_D236.mat
n=10^5;

%% 2.a
% Crude Monte Carlo
crude_montecarlo_means = zeros(length(lambda_vec),1);
crude_montecarlo_variances = zeros(length(lambda_vec),1);
crude_ci = zeros(length(lambda_vec),2);
crude_all_vals = zeros(length(lambda_vec),n);


for i = 1:length(lambda_vec)
    % Draw n X_i from weibull distribution
    wind = wblrnd(lambda_vec(i), k_vec(i), n, 1);
    % Calculate psi(X_i)
    power=P(wind);
    % Put all values for psi(X) in avector
    crude_all_vals(i,:) = power;
    % Take mean of all psi(X) and store in vector
    mean_power =  mean(power);
    crude_montecarlo_means(i) = mean_power;
    % Calculate variance and store
    var_power = var(power);
    crude_montecarlo_variances(i) = var_power;
    
    % Create CL (two sided) by using fact from  L2 p.6
    cl = 1.96*sqrt(var_power/n); 
    crude_ci(i,1) = mean_power + cl;
    crude_ci(i,2) = mean_power-cl;
end
figure(1)
boxplot(crude_all_vals')
title('Raw Monte Carlo')
ylabel('Power')
xlabel('Month')

% Truncated Monte Carlo
truncated_montecarlo_means = zeros(length(lambda_vec),1);
truncated_montecarlo_variances = zeros(length(lambda_vec),1);
truncated_ci = zeros(length(lambda_vec),2);
truncated_all_vals = zeros(length(lambda_vec),n);
for i =1:length(lambda_vec)
    % Weibull inv according to formula in 1.
    u = rand(n,1);
    [F_c_inv, c]  = truncate_wbl_dist(lambda_vec(i), k_vec(i),n,u);

    % New power 
    power = P((F_c_inv)).*(1-c)';
    truncated_all_vals(i,:) = power; %här multiplicerar vi med c
    mean_power = mean(power);
    truncated_montecarlo_means(i) = mean_power;
    
    var_power = var(power);
    truncated_montecarlo_variances(i) = var_power;
    cl = 1.96*sqrt(var_power/n); 
    truncated_ci(i,1) = mean_power + cl;
    truncated_ci(i,2) = mean_power-cl;
end 
figure(2)
boxplot(truncated_all_vals')
title('Truncated Monte Carlo')
ylabel('Power')
xlabel('Month')
%% 2.b Use wind V as control variate (CV = control variate)
% I rapporten måste vi vara noga med att förklara att korrelationen finns!
% Samt ta med lite bevis för Beta osv!

% här tar jag vår trunkerade varians och reducerar den!
CV_montecarlo_means = zeros(length(lambda_vec),1);
CV_montecarlo_variances = zeros(length(lambda_vec),1);
CV_ci = zeros(length(lambda_vec),2);
CV_all_vals = zeros(length(lambda_vec),n);
corrs = zeros(length(lambda_vec),1);
for i =1:length(lambda_vec)
    wind = wblrnd(lambda_vec(i), k_vec(i), n, 1);
    power = P(wind);
    mean_power = mean(power);
    
    old_var = var(power);
    % här kommer CV in!
    new_var = (1-corr(wind,power).^2)*old_var;
    corrs(i)  = corr(wind,power);
    CV_montecarlo_means(i) = mean_power;
    CV_montecarlo_variances(i) = new_var;
    cl = 1.96*sqrt(new_var/n);
    CV_ci(i,1) = mean_power + cl;
    CV_ci(i,2) = mean_power-cl;
end 
mean(corrs)
%% 2c Test with g as gamma dist
% IS = importance sampling
IS_montecarlo_means = zeros(length(lambda_vec),1);
IS_montecarlo_variances = zeros(length(lambda_vec),1);
IS_ci = zeros(length(lambda_vec),2);
IS_all_vals = zeros(length(lambda_vec),n);
for i =1:length(lambda_vec)
    % trunkerar gamma-fördelningen
    u = rand(n,1);
    [g_rnd, c] = truncate_gam_dist(9.7,1.25, n, u); 
    f = wblpdf(g_rnd, lambda_vec(i), k_vec(i));
    g = gampdf(g_rnd,9.7,1.25); 
    w=f./g;
    tau=P(g_rnd).*w'.*(1-c);
    IS_all_vals(i,:) = tau;
    mean_power=mean(tau);
    var_power=var(tau);
    
    IS_montecarlo_means(i) = mean_power;
    IS_montecarlo_variances(i) = var_power;
    cl = 1.96*sqrt(var_power/n);
    IS_ci(i,1) = mean_power + cl;
    IS_ci(i,2) = mean_power-cl;
end
figure(3)
boxplot(IS_all_vals')
title('Importance sampling')
ylabel('Power')
xlabel('Month')


%% 2.d
% Antithetic sampling ( = AS ) Obs vi måste använda truncated wind annars
% håller inte theorem!
AS_montecarlo_means = zeros(length(lambda_vec),1);
AS_montecarlo_variances = zeros(length(lambda_vec),1);
AS_ci = zeros(length(lambda_vec),2);
AS_all_vals = zeros(length(lambda_vec),n);
for i = 1:length(k_vec)
    
  
    u1 = rand(n,1);
    u2 = 1-u1;
    % truncate wind
    [F_c_inv, c1] = truncate_wbl_dist(lambda_vec(i), k_vec(i),n,u1);
    [F_c_inv2,c2] = truncate_wbl_dist(lambda_vec(i), k_vec(i),n,u2); % Vi använder inte truncated wind?
    
    V1 = P(F_c_inv).*(1-c1);
    V2 = P(F_c_inv2).*(1-c2);
   
    W = (V1+V2)/2;
    AS_all_vals(i,:) = W;
    tau_AS = mean(W);
    AS_montecarlo_means(i) = tau_AS;
    AS_montecarlo_variances(i)  = var(W);
    AS_ci(i,1) = tau_AS + norminv(0.975)*std(W)./sqrt(n/2); % Varför n/2?
    AS_ci(i,2)= tau_AS - norminv(0.975)*std(W)./sqrt(n/2);
end
figure(4)
boxplot(AS_all_vals')
title('Antithetic sampling')
ylabel('Power')
xlabel('Month')
%% 2.e

non_zero_prob = zeros(length(k_vec),1);
for i = 1:length(k_vec)
    wind = wblrnd(lambda_vec(i), k_vec(i), n, 1);
    power=P(wind);
    non_zero_prob(i) = nnz(power)/n;
end

figure(13)
plot(non_zero_prob, 'o')
title('Probability of positive power-output')
ylabel('Probability')
xlabel('Month')

%% 2.f

rho = 1.225;
d = 236;
coeff = 0.25*0.5*rho*pi*d^2;
ratio_means = zeros(length(lambda_vec),1);
ratio_variances = zeros(length(lambda_vec),1);
ratio_all_vals = zeros(length(lambda_vec),n);
for i=1:length(k_vec)
    E_p_tot = gamma(1+3/k_vec(i))*lambda_vec(i)^3;
    ratio = AS_all_vals(i,:) / (coeff*E_p_tot);
    ratio_all_vals(i,:) = ratio;
    ratio_variances(i) = var(AS_montecarlo_means(i))/(coeff*E_p_tot)^2;
    ratio_means(i) = mean(ratio);
    
end

figure(5)
boxplot(ratio_all_vals')
title('Ratio')
ylabel('Power ratio')
xlabel('Month')
%% 2.g

availability_vec = zeros(length(k_vec),1);
capacity_vec = zeros(length(k_vec),1);
for i=1:length(k_vec)
    k_temp = k_vec(i);
    lambda_temp = lambda_vec(i);
    wind_temp = wblrnd(lambda_temp, k_temp, n, 1);
    power_temp = P(wind_temp);

    max_capacity = 15*10^6;
    availability_vec(i) = nnz(power_temp)/n;
    capacity_vec(i) = crude_montecarlo_means(i) /max_capacity;
end

mean(capacity_vec)
mean(availability_vec)
figure(6) 
subplot(211)
plot(availability_vec,'o')
ylabel('Availability')
xlabel('Month')
subplot(212)
plot(capacity_vec,'ro')
ylabel('Capacity')
xlabel('Month')

%% Misc: plotta alla varianser och means
%kan vara trevligt att ha i rapporten

crude_var_mean = crude_montecarlo_variances';
truncated_var_mean = truncated_montecarlo_variances';
CV_var_mean = CV_montecarlo_variances';
IS_var_mean = IS_montecarlo_variances';
AS_var_mean = AS_montecarlo_variances';
figure(7)
plot([crude_var_mean truncated_var_mean CV_var_mean IS_var_mean AS_var_mean],'x')
title('Variance for different methods')
xticks([6 18 30 42 54])
xticklabels({'Raw' 'Trunc' 'CV' 'IS' 'AS'})
ylabel('Log-scale')
set(gca, 'YScale', 'log')


crude_mean = crude_montecarlo_means';
truncated_mean = truncated_montecarlo_means';
CV_mean = CV_montecarlo_means';
IS_mean = IS_montecarlo_means';
AS_mean = AS_montecarlo_means';
figure(8)
plot([crude_mean truncated_mean CV_mean IS_mean AS_mean],'x')
title('Means for different methods')
xticks([6 18 30 42 54])
xticklabels({'Raw' 'Trunc' 'CV' 'IS' 'AS'})
ylabel('Log-scale')
set(gca, 'YScale', 'log')


%% 3 a old
% to get expected value 
lambda = 10.05;
k = 1.95;

u = rand(n,1);
F_a = wblcdf(3,lambda,k);
F_b = wblcdf(30,lambda,k);
c = wblcdf(3,lambda,k) + (1- wblcdf(30,lambda,k));
% Weibull inv according to formula in 1.
x = wblinv(u.*(F_b-F_a)+F_a,lambda,k)'; 

f = wblpdf(x,lambda, k);
g = gampdf(x,9.70,1.25);
w =  f./ g;
mean_P1_P2 = mean(2*P(x) .*w');
var_P1_P2 = var(P(x).*w'); %samma för båda

%% 3b old
% covariance cov(X,Y) = E(XY)-E(X)E(Y)
mean_P1P2 = mean(P(x).*P(x).*w'); % Måste vi inte importance sample?
var_P1P2 = var(P(x).*P(x).*w');
cov_P1_P2 = mean_P1P2-(mean_P1_P2);

%% 3c old
% V(X+Y) = V(X) + V(Y) + 2C(X,Y)
% D(X+Y) = sqrt(V(X+Y))
var_P1_P2 = 2*var_P1_P2 + 2*cov_P1_P2;
std_P1_P2 = sqrt(var_P1_P2);

%% 3b,c new
% covariance cov(X,Y) = E(XY)-E(X)E(Y)
% Define wind and params
v1=0:0.5:30; 
v2=0:0.5:30;
lambda = 10.05;
k = 1.95;
V_1 = wblrnd(lambda, k, n, 1); % Generate separate winds
V_2 = wblrnd(lambda, k, n, 1);
alpha = 0.638;
p=3;
q=1.5;

wind_dist_pdf = @(v) wblpdf(v, 10.05, 1.95)';
wind_dist_cdf = @(v) wblcdf(v, 10.05, 1.95)';
f = @(v_1, v_2) wind_dist_pdf(v_1).*wind_dist_pdf(v_2).*(1+alpha*((1-wind_dist_cdf(v_1).^p).^(q-1)).*((1-wind_dist_cdf(v_2).^p).^(q-1)).*(wind_dist_cdf(v_1).^p.*(1+p*q)-1).*(wind_dist_cdf(v_2).^p.*(1+p*q)-1));
mu = [12, 12];
Sigma = [17, 3; 3, 17];
% Plot the distributions
[X,Y] = meshgrid(v1,v2);
px = reshape(P(X),length(v2),length(v1));
py = reshape(P(Y),length(v2),length(v1));
Z1=px.*py.*f(X,Y);
figure
surf(X,Y,Z1)
X = [X(:) Y(:)];
g = mvnpdf(X,mu,Sigma);
g = reshape(g,length(v2),length(v1));
Z2=g;
figure
surf(v1,v2,Z2);

% IS to get the estimates.
X1=mvnrnd(mu,Sigma,n);
eval=f(X1(:,1),X1(:,2))'; % Utvärdera f motsvarande
g_eval=mvnpdf(X1,mu,Sigma);
w=eval./g_eval;
%plot(w)
px=P(X1(:, 1)).*w;
py=P(X1(:, 2)).*w;

% Calculate covar, variance and std.
pxpy=px.*py;
cova = mean(pxpy) - mean(px).*mean(py);
cov_scalar=mean(cova);
variance = var(px)+var(py)+2*cov_scalar
std=sqrt(variance)

% Try with gamma dist - worse!
g = @(v1, v2) gampdf(v1, 9.7, 1.25) .* gampdf(v2, 9.7, 1.25);
g_rand1=gamrnd(9.7, 1.25,1, n)';
g_rand2=gamrnd(9.7, 1.25,1, n)';
eval=f(g_rand1,g_rand2)'; % Utvärdera f motsvarande
g_eval=g(g_rand1, g_rand2);
w=eval./g_eval;
%plot(w)
px=P(g_rand1).*w;
py=P(g_rand2).*w;
pxpy=px.*py.*w;
cova = mean(pxpy) - mean(px).*mean(py);
cov_scalar=mean(cova);
variance = var(px)+var(py)+2*cov_scalar;
std=sqrt(variance);

%% 3d part 1.
v1=0:0.5:35; % wind
v2=0:0.5:35;
wind_dist_pdf = @(v) wblpdf(v, 10.05, 1.95)'; % Params given
wind_dist_cdf = @(v) wblcdf(v, 10.05, 1.95)';
alpha = 0.638;
p=3;
q=1.5;
% Define multivariate density function:
f = @(v_1, v_2) wind_dist_pdf(v_1).*wind_dist_pdf(v_2).*(1+alpha*((1-wind_dist_cdf(v_1).^p).^(q-1)).*((1-wind_dist_cdf(v_2).^p).^(q-1)).*(wind_dist_cdf(v_1).^p.*(1+p*q)-1).*(wind_dist_cdf(v_2).^p.*(1+p*q)-1));
indicator1=@(v1, v2) P(v1)+P(v2)>15000000; % create indicator1

% To use IS we need a function g(x,y) to mimic phi(x,y)f(x,y)
% For the first indicator the objective function is defined on a finite
% grid. Therefore, we should be fine with sampling from a multivariate
% normal as before. But let's first plot that!

[X,Y] = meshgrid(v1,v2);
% Params for normal
mu = [12, 12];
Sigma = [35, 3; 3, 35];
X2 = [X(:) Y(:)];
g = mvnpdf(X2,mu,Sigma);
g = reshape(g,length(v2),length(v1));
ind=indicator1(X, Y);
ind=reshape(ind,length(v2),length(v1));
eval=f(X,Y).*ind;
figure
surf(v1,v2,eval); % Plot objective function times f

figure
surf(v1,v2,g) % Plot g

figure
surf(v1,v2,eval./g) % Plot ratio, bounded

% Use this and sample from g
X1=mvnrnd(mu,Sigma,n);
g_eval=mvnpdf(X1,mu,Sigma); % Evaluate all g
ind=indicator1(X1(:, 1), X1(:, 2)); % Evaluate all generated winds.
eval=f(X1(:,1),X1(:,2))'; % Evaluate f the same way

val=ind.*eval./g_eval; % Probability vector power>15MW.
probability=mean(val); % Expected value of above.
cl_over = 1.96*sqrt(var(val)/n); % 95% CL
total_prob=mean(eval./g_eval); % Total probability
%% 3 d part 2
% When evaluating the other probability, now the region is infinite. To be
% able to estimate the probability a new function g must be found that
% decays slower than f*phi. Here we try a 2d gamma function.

% Create new indicator as our objective function.
indicator2=@(v1, v2) P(v1)+P(v2)<15000000;

% Define 2d gamma function to replicate our indicator*f
const1 = 3.1;
const2 = 2.1;
g = @(v1, v2) gampdf(v1', const1, const2) * gampdf(v2, const1, const2);

% Plot to see if it seems resonable
[X,Y] = meshgrid(v1,v2);
g_wind=g(v1,v2);
figure
surf(v1,v2,g_wind)

figure
ind=indicator2(X, Y);
ind=reshape(ind,length(v2),length(v1));
eval=f(X,Y).*ind;
surf(v1, v2, eval)

figure
surf(v1,v2,eval./g_wind); % Good enough...

% Draw random variables from that distribution.
g_rand1=gamrnd(const1, const2, 1, n)';
g_rand2=gamrnd(const1, const2,1, n)';
% evaluate g, ind and f
g_eval=gampdf(g_rand1, const1, const2) .* gampdf(g_rand2, const1, const2); % evaluera alla värden i g
ind=indicator2(g_rand1, g_rand2); % Utvärdera alla genererade vindar
eval=f(g_rand1,g_rand2)'; % Utvärdera f motsvarande


val=ind.*eval./g_eval; % Probability vector P>15MW
probability=mean(val); % Expected value of above.
cl_under = 1.96*sqrt(var(val)/n); % 95% CL
total_prob=mean(eval./g_eval); % total probability

%% functions
function [t,c] = truncate_wbl_dist(lambda, k,n,u)
        F_a = wblcdf(3,lambda,k);
        F_b = wblcdf(30,lambda,k);
        c = wblcdf(3,lambda,k) + (1- wblcdf(30,lambda,k));
        t = wblinv(u.*(F_b-F_a)+F_a,lambda,k);
end

function [t,c] = truncate_gam_dist(lambda, k,n,u)
    F_a = gamcdf(3,lambda,k); 
    F_b = gamcdf(30,lambda,k);
    u = rand(1,n);
    t = gaminv(u.*(F_b-F_a)+F_a,lambda,k);
    c = gamcdf(3,lambda,k) + (1 - gamcdf(30,lambda, k));
    
end



