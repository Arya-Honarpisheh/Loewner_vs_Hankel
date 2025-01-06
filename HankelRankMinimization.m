function [h, p, v] = HankelRankMinimization(L, N, y, Tu_N, epsilon_t, rho)

%%%%%%%%%%%%   h = l11 + Tp*l12   %%%%%%%%%%%%%%
%%%%%%%%%%%%   v = l21 + Tp*l22   %%%%%%%%%%%%%%
%%%%%%%%%%%%         p = Qv       %%%%%%%%%%%%%%
disp("___________________________________________________________________")
disp("Hankel Method:")

l11 = impulse(minreal(L(1,1)), N-1); l12 = impulse(minreal(L(1,2)), N-1);
l21 = impulse(minreal(L(2,1)), N-1); l22 = impulse(minreal(L(2,2)), N-1);


Tl21 = toeplitz(l21, [l21(1), zeros(1,N-1)]);
Tl12 = toeplitz(l12, [l12(1), zeros(1,N-1)]);
Tl22 = toeplitz(l22, [l22(1), zeros(1,N-1)]);
p = sdpvar(N, 1, 'full', 'real');
Tp = toeplitz(p, [p(1), zeros(1, N-1)]);
% positive defeniteness of Qcond ensures us that Q is a stable system

% we need h to use it in the objective function
h = l11 + Tl12*p;
% Hankel matrix is nH*mH, [h(2); ... ; h(nH+1)] is the first column
% [h(nH+1), ... , h(N)] is the last row, and N = nH+mh
nH = floor(N/2); mH = N - nH;
Hro = hankel(h(2:nH+1)', h(nH+1:end));

%%%%%%%%% Reweighted Nuclear Norm Minimization %%%%%%%%%%%
Options = sdpsettings('solver', 'mosek', 'verbose',0, 'debug',0, 'savesolverinput',1);
R2i = diag(rho.^(-2*(0:N-1)));
% (RTvRi)(RTvRi)' - (RTpRi)(RTpRi)' >= 0
Qcond = [Tl21*R2i*Tl21' + Tl21*R2i*Tl22'*Tp' + Tl22*Tp*R2i*Tl21',   Tp; Tp', inv(R2i*eye(N) - Tl22*R2i*Tl22')];
Qcond = (Qcond+Qcond')/2;
ConstraintsHankel = [Qcond >= 0, norm(y-Tu_N*h, 'inf') <= epsilon_t];

counter = 0;
counter_max = 5;
delta = 0.001;
W1 = eye(nH); W2 = eye(mH);
Y = eye(nH); Z = eye(mH);
while( counter < counter_max )
    Objective = norm(W1*Hro*W2, '*');
    optimize(ConstraintsHankel, Objective, Options);
    counter = counter + 1;
    [U, S, V] = svd(W1*value(Hro)*W2, "econ");
    Y = (Y + delta*eye(nH))^(1/2) * U * S * U' * (Y + delta*eye(nH))^(1/2);
    Z = (Z + delta*eye(mH))^(1/2) * V * S * V' * (Z + delta*eye(mH))^(1/2);
    W1 = (Y + delta*eye(nH))^(-1/2);
    W2 = (Z + delta*eye(mH))^(-1/2);
end

fprintf("The smallest eigen value of Qcond is %f >= 0 \n", min(eig(value(Qcond))))
fprintf("Consistensy infimum norm is %f <= %f \n", norm(y-Tu_N*value(h), 'inf'), epsilon_t)
disp("___________________________________________________________________")

h = value(h);
p = value(p);
v = l21 + value(Tp)*l22;

end