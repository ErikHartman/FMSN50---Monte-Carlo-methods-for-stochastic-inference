%% 3

% using a naive sequential importance sampling
N = 100;
d = 2;
k_max = 2;
X = zeros(k_max,d, N);
w = zeros(k_max, N);
X_SA = zeros(1,N);
for n = 1:N
    for k = 1:k_max
        curr_pos = X(k,:,n);
        possible_new_pos = get_positions(curr_pos);
        random_index = randi(size(possible_new_pos,1));
        new_pos = possible_new_pos(random_index,:);
        X(k+1,:,n) = new_pos;
    end
    X_SA(n) =  length(unique(X(:,:,n),'rows')) == k_max+1;
end
ratio = nnz(X_SA) / N


%% 4

N_max = 10;
d = 2;
k_max = 100;
estimated_c = zeros(N_max,1);
for N = 2:N_max
    X = zeros(k_max,d, N);
    c = zeros(N,1);
    w = zeros(k_max, N);
    w(1,:) = 1;
    X_SA = zeros(1,N);
    for k = 1:k_max
        for n = 1:N
            history = X(:,:,n);
            curr_pos = X(k,:,n);
            possible_new_pos = get_positions(curr_pos);
            possible_new_pos = get_valid_positions(history, possible_new_pos);
            if isempty(possible_new_pos)
                X(k+1,:,n) = curr_pos;
                w(k+1,n) = 0;
            else
                if(size(possible_new_pos,1) > 3)
                     w(k+1,n) = 1;
                else
                    w(k+1,n) = (size(possible_new_pos,1)/3) * w(k,n);
                end
                random_index = randi(size(possible_new_pos,1));
                new_pos = possible_new_pos(random_index,:);
                X(k+1,:,n) = new_pos;
            end
        end
        w = normalize_weights(w);
    end
    for n = 1:N
        [~,l,~] =  unique(X(:,:,n),'rows');
        c(n) = max(l);
        weight_at_index = w(max(l),n);
        c(n) = c(n) .* weight_at_index;
        
    end
    E_c = sum(c);
    estimated_c(N) = E_c;
    
end
SIS_est = mean(estimated_c)
figure(1)
plot(estimated_c, 'o')
title('SIS - estimated c')
xlabel('N')
ylabel('c')

%% 5
clear
N_max = 100;
d = 2;
k_max = 50;
estimated_c = zeros(N_max,1);

for N = 20:N_max
    X = zeros(k_max,d, N);
    c = zeros(N,1);
    w = zeros(k_max, N);
    w(1,:) = 1;
    X_SA = zeros(1,N);
    for k = 1:k_max
        for n = 1:N
            history = X(:,:,n);
            curr_pos = X(k,:,n);
            possible_new_pos = get_positions(curr_pos);
            possible_new_pos = get_valid_positions(history, possible_new_pos);
            if isempty(possible_new_pos)
                X(k+1,:,n) = curr_pos;
                w(k+1,n) = 0;
            else
                if(size(possible_new_pos,1) > 3)
                     w(k+1,n) = 1;
                else
                    w(k+1,n) = (size(possible_new_pos,1)/3) * w(k,n);
                end
                random_index = randi(size(possible_new_pos,1));
                new_pos = possible_new_pos(random_index,:);
                X(k+1,:,n) = new_pos;
            end
        end

        w(isnan(w))=0;
        if(sum(w(k,:) ~= 0))
            [X,w] = resample(X,w,N,k);
            w = normalize_weights(w);
        end
    end
    for n = 1:N
        [~,l,~] =  unique(X(:,:,n),'rows');
        c(n) = max(l);
        weight_at_index = w(max(l),n);
        c(n) = c(n) .* weight_at_index;       
    end
    E_c = sum(c);
    estimated_c(N) = E_c;
end

SIRS_est = mean(estimated_c)
figure(2)
plot(estimated_c, 'o')
title('SIRS - estimated c')
xlabel('N')
ylabel('c')
%% functions

function [X,w] = resample(X,w,N,k)
    ind = randsample(N,N,true,w(k,:));
    X = X(:,:, ind);
    w(k,ind) = 1; 
end

function w = normalize_weights(weights)
    w = weights ./ sum(weights,2);
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


