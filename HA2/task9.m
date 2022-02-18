

tries = 10;

N = 100;
d = 100;
k_max = 20;
A = zeros(tries,1);
mu = zeros(tries,1);
gamma = zeros(tries,1);

for t = 1:tries
    t
    X = zeros(k_max,d, N);
    SA = zeros(k_max,1);
    w = zeros(k_max,N);
    w(1, :) = 1;
    for k = 1:k_max
        for n = 1:N
            curr_pos = X(k,:,n);
            history = X(:,:,n);
            possible_new_pos = get_valid_positions2(X,k,n);
            if isempty(possible_new_pos)
                X(k+1,:,n) = curr_pos;
                w(k+1,n) = 0;
            else
                random_index = randi(size(possible_new_pos,1));
                X(k+1,:,n) = possible_new_pos(random_index,:);
                w(k+1,n) = w(k,n) * (length(possible_new_pos)); 
            end
        end
    end
    cn = mean(w(2:end,:),2);
    Y = log(cn) + log(1:k_max)';
    X_reg = [ones(k_max,1) (1:k_max)' log(1:k_max)'];
    beta = X_reg\Y;
    expbeta = exp(beta)';

    A2_reg = expbeta(1);
    mu2_reg = expbeta(2);
    gamma2_reg = beta(3);

    A(t) = A2_reg;
    mu(t) = mu2_reg;
    gamma(t) = gamma2_reg;
end
mean(mu)
Graham = 2*d-1-1/(2*d)-3/(2*d)^2-16/(2*d)^3

function possible_nodes = get_valid_positions2(X,step,n)
    dim = size(X);
    dim = dim(2);
    all_directions = cat(2,eye(dim),-eye(dim));
    current_node = X(step,:,n)';
    all_new_nodes = (current_node + all_directions)';
    possible_nodes = [];
    for i=1:2*dim
        if ~any(ismember(X(:,:,n), all_new_nodes(i,:), 'rows')); % Directions not in history
            possible_nodes(end+1,:)=all_new_nodes(i,:); % If not add to possible nodes
        end
    end
    
end