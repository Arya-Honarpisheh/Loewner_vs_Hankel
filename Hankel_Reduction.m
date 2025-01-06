

function G = Hankel_Reduction(h, tol)
% Ho-Kalman algorithm

N = length(h);
nH = floor(N/2); mH = N - nH;
H = hankel(h(2:nH+1)', h(nH+1:end));

Hl = H(:,1:mH-1); Hr = H(:,2:end);
[U, S, V] = svd(Hl);
red_ord = find(cumsum(diag(S)) / sum(diag(S)) > tol, 1);
S = S(1:red_ord, 1:red_ord);
Kl = U(:,1:red_ord)*S^(1/2); % Kl is the observability matrix
Kr = S^(1/2)*V(:,1:red_ord)'; % Kr is the contolability matrix
A = pinv(Kl) * Hr * pinv(Kr);
B = Kr(:,1);
C = Kl(1,:);
D = h(1);
G = minreal(ss(A, B, C, D, 1), 1e-3);
G = tf(G);

end