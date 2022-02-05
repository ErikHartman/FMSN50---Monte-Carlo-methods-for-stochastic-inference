

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

for i =1:length(lambda_vec)
    u = rand(n,1);
    [F_c_inv, c] = truncate_wbl_dist(lambda_vec(i), k_vec(i), n,u);
    power = P((F_c_inv)).*(1-c);
    mean_power = mean(power);
    
    old_var = var(power);
    % här kommer CV in!
    new_var = (1-corr(F_c_inv,power).^2)*old_var;
    
    CV_montecarlo_means(i) = mean_power;
    CV_montecarlo_variances(i) = new_var;
    cl = 1.96*sqrt(new_var/n);
    CV_ci(i,1) = mean_power + cl;
    CV_ci(i,2) = mean_power-cl;
end 

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
nnz(power)/n;

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
    capacity_vec(i) = AS_montecarlo_means(i) /max_capacity; % not our best estimate
end

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

%% 3 a
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

%% 3b
% covariance cov(X,Y) = E(XY)-E(X)E(Y)
mean_P1P2 = mean(P(x).*P(x).*w'); % Måste vi inte importance sample?
var_P1P2 = var(P(x).*P(x).*w');
cov_P1_P2 = mean_P1P2-(mean_P1_P2);

%% 3c
% V(X+Y) = V(X) + V(Y) + 2C(X,Y)
% D(X+Y) = sqrt(V(X+Y))
var_P1_P2 = 2*var_P1_P2 + 2*cov_P1_P2;
std_P1_P2 = sqrt(var_P1_P2);

%% 3d
% Vi börjar med att kolla hur phi(x)*f(x,y) ser ut för att hitta en lämplig
% multivariate g(x,y)
v1=0:0.5:35; % vind i vanlig ordning
v2=0:0.5:35;
wind_dist_pdf = @(v) wblpdf(v, 10.05, 1.95)'; % Parametrar givna i uppgiften
wind_dist_cdf = @(v) wblcdf(v, 10.05, 1.95)';
alpha = 0.638;
p=3;
q=1.5;

%Multivariate density function:
f = @(v_1, v_2) wind_dist_pdf(v_1).*wind_dist_pdf(v_2).*(1+alpha*((1-wind_dist_cdf(v_1).^p).^(q-1)).*((1-wind_dist_cdf(v_2).^p).^(q-1)).*(wind_dist_cdf(v_1).^p.*(1+p*q)-1).*(wind_dist_cdf(v_2).^p.*(1+p*q)-1));

% Kör 3D plottar, här plottas phi(x)*phi(y)*f(x,y) (är ej helt säker om det
% räcker med en av dom...
[X,Y] = meshgrid(v1,v2);
px = reshape(P(X),length(v2),length(v1));
py = reshape(P(Y),length(v2),length(v1));
Z1=f(X,Y).*px.*py;
figure
surf(X,Y,Z1); % Ser ut typ som en multivariant normalfördelning kan passa ok!
xlabel('v1')
ylabel('v2')

% Bygg multivariant normal, testade mig fram med parametrarna. Verkar som
% ju mer off desto mer påverkas total sannolikhet, typ rimligt?
mu = [12, 12];
Sigma = [16, 3; 3, 16];

X = [X(:) Y(:)];
g = mvnpdf(X,mu,Sigma);
g = reshape(g,length(v2),length(v1));
Z2=g;
figure
surf(v1,v2,Z2);

% Generera slumptal från den fördelningen.
X1=mvnrnd(mu,Sigma,n);
g_eval=mvnpdf(X1,mu,Sigma); % evaluera alla värden i g

indicator=@(v1, v2) P(v1)+P(v2)>15000000; % skapa indikator, ny phi(x)
ind=indicator(X1(:, 1), X1(:, 2)); % Utvärdera alla genererade vindar
f_eval=f(X1(:,1),X1(:,2))'; % Utvärdera f motsvarande
val=ind.*f_eval./g_eval; % Sannolikhetsvektor att P>15MW
probability=mean(val) % Sannolikheten för ovanstående, runt 0.5
total_prob=mean(f_eval./g_eval) % total sannolikhet (utan indikator), runt 1



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



