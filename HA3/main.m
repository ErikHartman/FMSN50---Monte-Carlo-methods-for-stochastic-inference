%% Part 1
load('coal_mine_disasters.mat')

% Plot the data
histogram (T)

% Define initial gamma(2,psi) to sample our theta, start by setting psi=1.
gam_theta = @(psi,n) gamrnd(2,psi, n,1)