%% Task 2.
load('atlantic.txt')
n = length(atlantic);
[beta,mu] = est_gumbel(atlantic);

nsamples = 1000;
% Bootstrap parameters
bs_mus = zeros(nsamples,1);
bs_betas = zeros(nsamples,1);
bs_retwave = zeros(nsamples,1);

T = 3*14*100; % given
big_wave = Finv(1 - 1/T, mu, beta);
% Bootstrap
for i = 1:nsamples
    u = rand(n,1);
    fake_atlantic = Finv(u,mu,beta);
    [bs_beta, bs_mu] = est_gumbel(fake_atlantic);
    bs_retwave(i) = Finv(1 - 1/T, bs_mu, bs_beta);
    bs_mus(i) = bs_mu;
    bs_betas(i) = bs_beta;
end


% Calculating deltas
delta_wave =  big_wave - bs_retwave;
delta_mu = mu - bs_mus;
delta_beta = beta - bs_betas;
delta_mu = sort(delta_mu);
delta_beta = sort(delta_beta);
delta_wave = sort(delta_wave);

% Calculating confidence intervals
alpha = 0.05; % confidence 
conf_mu =  mu - [delta_mu(ceil((1 - alpha/2)*nsamples)), delta_mu(ceil(alpha*nsamples/2))]
conf_beta =  beta - [delta_beta(ceil((1 - alpha/2)*nsamples)), delta_beta(ceil(alpha*nsamples/2))]
conf_wave = big_wave + delta_wave(ceil((1 - alpha)*nsamples))


function k = Finv(u,mu,beta)
    k = mu-beta*log(-log(u));
end

