function HC = Eigen_neg_to_pos(Q)
N = size(Q,1);
HC = Q;
[U, E] = eig(Q,'vector');
v = E(E < 0);
m = length(v); % number of negative values
if m > 0
    S = sum(v);
    W = (S*S*100)+1;
    P = abs(E(m)); %smallest absolute value of negative eignvalue
    for i = 1:m
        C = E(i);
        E(i) = P * (S-C)*(S-C)/W;
    end
    HC = U * diag(E) * U';
end
%EOF
end