function [z, w, h, p, v, dcg] = LoewnerRankMinimization(L, N, y, Tu_N, epsilon_t, rho)

%%%%%%%%%%%%       w = L11(z) + L12(p)   %%%%%%%%%%%%%%
%%%%%%%%%%%%       v = L21(z) + L22(p)   %%%%%%%%%%%%%%
%%%%%%%%%%%%       p = Qv                %%%%%%%%%%%%%%
disp("___________________________________________________________________")
disp("Loewner Method:")

L11 = minreal(L(1,1)); L12 = minreal(L(1,2));
L21 = minreal(L(2,1)); L22 = minreal(L(2,2));

n = floor(N/2);
% theta = logspace(log10(pi/(n+1)), log10(pi-pi/(n+1)), n)';
theta = linspace(pi/(n+1), pi-pi/(n+1), n)';
za = exp(theta*1i);
zb = conj(za);
z = [za; zb];

FTM = zeros(N,N);
for i=1:N
    for j=0:N-1
        FTM(i,j+1) = z(i)^(-j);
    end
end

Fl11 = zeros(N, 1);
Fl12 = zeros(N, 1);
Fl21 = zeros(N, 1);
Fl22 = zeros(N, 1);
for k=1:N
    Fl11(k) = evalfr(L11, z(k));
    Fl12(k) = evalfr(L12, z(k));
    Fl21(k) = evalfr(L21, z(k));
    Fl22(k) = evalfr(L22, z(k));
end

p = sdpvar(N, 1, 'full', 'complex');
Z22 = diag(Fl22);
Qcond = [rho^-2*(Fl21'*Fl21) + rho^-2*Fl21'*Z22*p + rho^-2*p'*Z22'*Fl21, p'; p, inv(eye(N) - rho^-1*(Z22'*Z22)*rho^-1)];
% eig(eye(N) - rho*(Z22'*Z22)*rho)
Qcond = (Qcond+Qcond')/2;
% positive defeniteness of Qcond ensures us that Q is a stable system
% with Hinf norm less than one

w = Fl11 + Fl12.*p;
wa = w(1:n); wb = w(n+1:end);
Lro = sdpvar(n, n, 'full');
for i = 1:n
    for j = 1:n
        Lro(i,j) = (wb(i) - wa(j)) / (zb(i) - za(j));
    end
end

% to check the consistency with the data
dcg = sdpvar(1, 1, 'full', 'real');
h = FTM\w + dcg;

%%%%%%%%% Reweighted Nuclear Norm Minimization %%%%%%%%%%%
Options = sdpsettings('solver', 'mosek', 'verbose',0, 'debug',1);
ConstraintsLoewner = [Qcond >= 0 , norm(y-Tu_N*h, 'inf') <= epsilon_t];
% ConstraintsLoewner = [norm(y-Tu_N*h, 'inf') <= epsilon_t];

counter = 0;
counter_max = 5;
delta = 0.01;
W1 = eye(n); W2 = eye(n);
Y = eye(n); Z = eye(n);
while( counter < counter_max )
    Objective = norm(W1*Lro*W2, '*');
    optimize(ConstraintsLoewner, Objective, Options);
    counter = counter + 1;
    [U, S, V] = svd(W1*value(Lro)*W2, "econ");
    Y = (Y + delta*eye(n))^(1/2) * U * S * U' * (Y + delta*eye(n))^(1/2);
    Z = (Z + delta*eye(n))^(1/2) * V * S * V' * (Z + delta*eye(n))^(1/2);
    W1 = (Y + delta*eye(n))^(-1/2);
    W2 = (Z + delta*eye(n))^(-1/2);
end

fprintf("The smallest eigen value of Qcond is %f >= 0 \n", min(eig(value(Qcond))))
fprintf("Consistensy infimum norm is %f <= %f \n", norm(y-Tu_N*value(h), 'inf'), epsilon_t)
disp("___________________________________________________________________")

z = za;
w = value(wa);
p = value(p);
h = real(value(h));
v = Fl21 + Z22*p;
dcg = sum(h);

end