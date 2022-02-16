%% 3

% using a naive sequential importance sampling
N = 1000;
d = 2;
k_max = 11;
X = zeros(k_max,d, N);
SA = zeros(k_max,1);
w = zeros(k_max,N);
w(1, :) = 1;
for k = 1:k_max
    for n = 1:N
        curr_pos = X(k,:,n);
        history = X(:,:,n);
        possible_new_pos = get_positions(curr_pos);
        free_pos = get_valid_positions(history, possible_new_pos);
        z = 1;
        if(length(unique(X(1:k,:,n),'row')) < k)
            z = 0;
        end
        % Set weight
        w(k+1,n) = w(k,n) * (z/(1/(length(free_pos)))); 
        random_index = randi(size(possible_new_pos,1));
        new_pos = possible_new_pos(random_index,:);
        
        X(k+1,:,n) = new_pos;
    end
    SA(k) = get_number_SA(X, k);
end
cn = mean(w(2:end,:),2);
SA/N
figure;
plot(cn)
xlabel('Walk length')
ylabel('Approximate number of SAW') 

%% 4 new

% using a naive sequential importance sampling
N = 1000;
d = 2;
k_max = 11;
X = zeros(k_max,d, N);
SA = zeros(k_max,1);
w = zeros(k_max,N);
w(1, :) = 1;
for k = 1:k_max
    for n = 1:N
        curr_pos = X(k,:,n);
        history = X(:,:,n);
        possible_new_pos = get_positions(curr_pos);
        possible_new_pos = get_valid_positions(history, possible_new_pos);
        if isempty(possible_new_pos)
            X(k+1,:,n) = curr_pos;
            w(k+1,n) = 0;
        else
            random_index = randi(size(possible_new_pos,1));
            X(k+1,:,n) = possible_new_pos(random_index,:);
            w(k+1,n) = w(k,n) * (1/(1/(length(possible_new_pos)))); 
        end
    end
    SA(k) = get_number_SA(X, k);
end
cn = mean(w(2:end,:),2);
SA/N
figure;
plot(cn)
xlabel('Walk length')
ylabel('Approximate number of SAW') 

%% 5 new

N = 10;
d = 2;
k_max = 5;
X = zeros(k_max,d, N);
SA = zeros(k_max,1);
w = zeros(k_max,N);
w(1, :) = 1;
for k = 1:k_max
    for n = 1:N
        curr_pos = X(k,:,n);
        history = X(:,:,n);
        possible_new_pos = get_positions(curr_pos);
        possible_new_pos = get_valid_positions(history, possible_new_pos);
        if isempty(possible_new_pos)
            X(k+1,:,n) = curr_pos;
            w(k+1,n) = 0;
        else
            random_index = randi(size(possible_new_pos,1));
            X(k+1,:,n) = possible_new_pos(random_index,:);
            w(k+1,n) = w(k,n) * (1/(1/(length(possible_new_pos)))); 
        end
    end
    [X,w] = resample(X,w,N,k); %borde
    SA(k) = get_number_SA(X, k); 
end
cn = cumprod(mean(w(2:end,:),2));
SA/N
figure;
plot(cn)
xlabel('Walk length')
ylabel('Approximate number of SAW') 


%% 6

Y = log(cn) + log(1:k_max)';
X_reg = [ones(k_max,1) (1:k_max)' log(1:k_max)'];

beta = X_reg\Y;
expbeta = exp(beta)';

A2_reg = expbeta(1);
mu2_reg = expbeta(2);
gamma2_reg = beta(3);
%% Task 9, general d
% using a naive sequential importance sampling
N = 1;
d = 3;
k_max = 11;
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
            w(k+1,n) = w(k,n) * (1/(1/(length(possible_new_pos)))); 
        end
    end
    SA(k) = get_number_SA(X, k);
end
cn = mean(w(2:end,:),2);
figure;
plot(cn)
%plot3(X(:,1,1), X(:,2,1),X(:,3,1), 'b') % Plot random walk in 3d
xlabel('Walk length')
ylabel('Approximate number of SAW') 
Y = log(cn) + log(1:k_max)';
X_reg = [ones(k_max,1) (1:k_max)' log(1:k_max)'];

beta = X_reg\Y;
expbeta = exp(beta)';

A2_reg = expbeta(1);
mu2_reg = expbeta(2);
gamma2_reg = beta(3);

%% Task 10
N = 100;
n = 50;
tau = zeros(1,n); % vector of filter means
w = zeros(N,1);
p = @(x,y) unifpdf(y, x/2, x); % observation density, for weights
part = unifrnd(1/5,3/5,N,1); % initialization
w = p(part, Y(1)); % weighting
tau(1) = sum(part.*w)/sum(w); % estimation
ind = randsample(N,N,true,w); % selection
part = part(ind);
for k = 1:n % main loop
    B = unifrnd(1/2,3,N,1); 
    part = B.*part.*(1-part) % mutation
    w = p(part, Y(k+1)) % weighting
    tau(k + 1) = sum(part.*w)/sum(w) % estimation
    ind = randsample(N,N,true,w) % selection
    part = part(ind)
    [xx,I]=sort(part); % sort data
    cw=cumsum(w(I))/sum(w); % cumulative normalized weightsum
    % for sorted data
    Ilower=find(cw>=0.025,1); % index for lower 2.5% quantile
    Iupper=find(cw>=0.975,1); % index upper 2.5% quantile
    taulower(k+1)=xx(Ilower); % lower 2.5% quantile
    tauupper(k+1)=xx(Iupper); % upper 2.5% quantile
end

%%
load('population.mat')
hold on
plot(tau, 'b')

plot(X, 'black')
plot(taulower, 'r')
plot(tauupper, 'r')
%plot(Y)
legend("Tau","X", "taulower","tauupper")

hold off

%% functions
function [X,w] = resample(X,w,N,k)
    ind = randsample(N,N,true,w(k,:));
    X = X(:,:, ind);
    w(k+1,ind) = 1; 
end

function possible_new_positions = get_positions(curr_pos)
    possible_new_positions = zeros(4,2);
    possible_new_positions(1,:) = [curr_pos(1)-1 curr_pos(2)];
    possible_new_positions(2,:) = [curr_pos(1)+1 curr_pos(2)];
    possible_new_positions(3,:) = [curr_pos(1) curr_pos(2)-1];
    possible_new_positions(4,:) = [curr_pos(1) curr_pos(2)+1]; 
end

function pos = get_valid_positions(old_history, possible_new_positions)
    pos = [];
    k = 1; 
    for i =1:4
        if ~any(ismember(old_history, possible_new_positions(i,:), 'rows'))
            pos(k,:) = possible_new_positions(i,:);
            k = k+1;
        end
    end
end

function n = get_number_SA(X, k_max)
    n = 0;
    for i=1:length(X(1,1,:))
        current_pop = X(:,:,i);
        is_unique = length(unique(current_pop,'rows')) == k_max+1;
        if is_unique
            n = n+1;
        end
    end

end

% Function take previous positions and returns a random direction of the
% allowed ones in any dimension
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



