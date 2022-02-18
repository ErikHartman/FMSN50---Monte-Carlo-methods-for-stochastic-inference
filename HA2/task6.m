
A = zeros(10,1);
mu = zeros(10,1);
gamma = zeros(10,1);


N = 1000;
k_max = 50 + 1; 
d = 2; 
X = zeros(k_max,d,N);
w = zeros(k_max,N);
tries = 10;
estimates = zeros(k_max-1,tries);
w(1,:) = 1;

for t = 1:tries
    t
    for k = 2:k_max
       wd = cumsum(w(k - 1,:)./sum(w(k - 1,:)));
       for n = 1:N
            parent_index = find((wd > rand) == 1);
            parent_index = parent_index(1);
            X(1:k - 1,:,n) = X(1:k - 1,:,parent_index);         
            [new, nPossible] = improved_possible_dir(X(1:k - 1,:,n));
            X(k,:,n) = X(k - 1,:,n) + new;
            z = 1;
            if(length(unique(X(1:k,:,n),'row')) < k)
                z = 0;
            end
            w(k,n) = (z/(1/(nPossible))); 
       end 
    end
    cn = cumprod(mean(w(2:end,:),2));
    Y = log(cn) + log(1:k_max-1)';
    X_reg = [ones(k_max-1,1) (1:k_max-1)' log(1:k_max-1)'];

    beta = X_reg\Y;
    expbeta = exp(beta)';

    A2_reg = expbeta(1);
    mu2_reg = expbeta(2);
    gamma2_reg = beta(3);

    A(t) = A2_reg;
    mu(t) = mu2_reg;
    gamma(t) = gamma2_reg;
end


A
mu
gamma


%% functions
function [newDir, nPossible] = improved_possible_dir(X)
dims = size(X,2);
% Check all directions
oneHotDirs = [eye(dims); -eye(dims)];
possibleDirs = zeros(2*dims, dims);
for i = 1:2*dims
    % If direction is free, add direction to possibleDirs
    if isempty(intersect(X, X(end,:) + oneHotDirs(i,:), 'rows'))
        possibleDirs(i,:) = oneHotDirs(i,:);
    end
end
% Remove the zero rows
possibleDirs = possibleDirs(~all(possibleDirs == 0,2), :);

% Count how many of them are free
nPossible = size(possibleDirs,1);
% Choose a free direction at random

if nPossible == 0
   newDir = zeros(1,dims); 
elseif nPossible == 1
   newDir = possibleDirs;
else
   newDir = possibleDirs(randi(nPossible), :);
end
end

