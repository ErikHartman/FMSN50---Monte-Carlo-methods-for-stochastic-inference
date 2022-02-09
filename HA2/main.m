%% 3

% using a naive sequential importance sampling


N = 20000;
d = 2;
k_max = 10;
X = zeros(k_max,d, N);
w = zeros(k_max, N);
X_SA = zeros(1,N);

old_direction = [-1 -1];
for n = 1:N
    x = 0;
    y = 0;
    for k = 2:k_max+1
        direction = round(rand(1,d));
        if k > 3
            while (direction == [-old_direction(1), old_direction(2)] )
                direction = round(rand(1,d));
            end
        end
        old_direction = direction;
        if direction(1) == 1 && direction(2) == 0
            y = y-1;
            X(k,1,n) = y;
        end
        if direction(1) == 0 && direction(2) == 0
            x = x-1;
            X(k,2,n) = x;
        end
        if direction(1) == 1 && direction(2) == 1
            y = y+1;
            X(k,1,n) = y;
        end
        if direction(1) == 0 && direction(2) == 1
            x = x+1;
            X(k,2,n) = x;
        end
    end
    X_SA(n) =  length(unique(X(:,:,n),'rows')) == k_max+1;
end
ratio = nnz(X_SA) / N


%% 4

N = 1000;
d = 2;
k_max = 500;
c = zeros(N,1);
X = zeros(k_max,d, N);
w = zeros(k_max, N);
X_SA = zeros(1,N);
w(1,:) = 1;

for n = 1:N
    for k = 1:k_max
        old_history = X(:,:,n);
        curr_pos = X(k,:,n);
        possible_new_pos = get_positions(old_history,curr_pos);
        if isempty(possible_new_pos)
            X(k+1,:,n) = curr_pos;
            w(k+1,n) = 0;
            c(n) = k;
            break
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
end

hist(c)
% plot(X(:,1,1),X(:,2,1))
% xlim([-10,10])
% ylim([-10,10])


function pos = get_positions(old_history, curr_pos)
    possible_new_positions = zeros(4,2);
    possible_new_positions(1,:) = [curr_pos(1)-1 curr_pos(2)];
    possible_new_positions(2,:) = [curr_pos(1)+1 curr_pos(2)];
    possible_new_positions(3,:) = [curr_pos(1) curr_pos(2)-1];
    possible_new_positions(4,:) = [curr_pos(1) curr_pos(2)+1];
    
    pos = [];
    k = 1; 
    
    for i =1:4
        if ~any(ismember(old_history, possible_new_positions(i,:), 'rows'))
            pos(k,:) = possible_new_positions(i,:);
            k = k+1;
        end
    end
end