
function [G, L, L_sft] = Loewner_Reduction(wa, wb, za, zb, dcg, tol)
% Similar to Ho-Kalman algorithm

n = length(za); N = 2*n;
z = [za; zb];
FTM = zeros(N,N);
for i=1:N
    for j=0:N-1
        FTM(i,j+1) = z(i)^(-j);
    end
end

L = zeros(n,n);
L_sft = zeros(n,n);
for i = 1:n
    for j = 1:n
        L(i,j) = (wb(i) - wa(j)) / (zb(i) - za(j));
        L_sft(i,j) = (zb(i)*wb(i) - za(j)*wa(j)) / (zb(i) - za(j));
    end
end
red_ord = find(cumsum(diag(svd(L))) / sum(diag(svd(L))) > tol, 1);

% [U, S, V] = svd(-L, "econ");
% r = red_ord;
% Sigma = S(1:r,1:r);
% T = [inv(Sigma), zeros(r,n-r); zeros(n-r,r), eye(n-r,n-r)];
% Anew = T*U'*(-L_sft)*V;
% A11 = Anew(1:r,1:r);
% A22 = Anew(r+1:end,r+1:end);
% A21 = Anew(r+1:end,1:r);
% A12 = Anew(1:r,r+1:end);
% Bnew = T*U'*wb;
% B1 = Bnew(1:r);
% B2 = Bnew(r+1:end);
% Cnew = wb'*V;
% C1 = Cnew(1:r);
% C2 = Cnew(r+1:end);
% A = A11 - A12*pinv(A22)*A21;
% B = B1 - A12*pinv(A22)*B2;
% C = C1 - C2*pinv(A22)*A21;
% D = -C2*pinv(A22)*B2;

% H = dlyap(L_sft, wb*wb'*1e-10, [], L);
% L_s*H*L_s - L*H*L = -wb*wb'
% svd(H)

% Gid = ss(A,B,C,D,-1);
yek = ones(n,1);
% implementing bisection to find d such that the dcgain(Gid) = d
low = -10; high = 10;
while abs(high - low) > 0.0001
    d = (low + high) / 2;
    Gid = dss(-L_sft+yek*d*yek', wb-yek*d, transpose(wa) - d*yek', d, -L, -1);
    Gid = tf(Gid);
    Gid = tf(real(Gid.num{1}),real(Gid.den{1}),-1);
    Gid = ss(Gid);
    Gid = minreal(stabsep(Gid));
    if dcgain(Gid) < dcg
        low = d;
    else
        high = d;
    end
end
% This system is consistent with the data set (z,w)
% Since wb=conj(wa) and L, Ls are Hermitian Gid has a representation
% as a state space model with real parameteres i.e. all the poles of the
% system are either real or complex conjugate
% red_ord = 4;
if length(Gid.A) > red_ord
    G = balred(minreal(Gid), red_ord);
else
    G = Gid;
end

end