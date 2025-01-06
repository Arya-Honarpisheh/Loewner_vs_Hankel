function L = CminCposApi2LFT(Cmin, Cpos, Api, rho, K, cnjgt)

N = length(Cmin);
Lambda = dlyap(Api', Cmin'*Cmin - Cpos'*Cpos);
% fprintf("The minimum eigen value of the Pick matrix is %f >= 0 \n", min(eig(Lambda)));
% disp("_________________________________________________________________")

% Now we make T
A_T = Api;
C_T = [Cpos; Cmin]*(Api - eye(N));
B_T = ((Api'-eye(N))*Lambda)\[-Cpos' Cmin'];
D_T = eye(2) + [Cpos; Cmin] * B_T;

% Now we convert T into the L such that LFT(L,Q) parmetrizes all the interpolants
A_L = A_T - B_T(:,2)*pinv(D_T(2,2))*C_T(2,:);
B_L = [B_T(:,2)*pinv(D_T(2,2)), B_T(:,1)-B_T(:,2)*pinv(D_T(2,2))*D_T(2,1)];
C_L = [(C_T(1,:) - D_T(1,2)*pinv(D_T(2,2))*C_T(2,:)); -pinv(D_T(2,2))*C_T(2,:)];
D_L = [D_T(1,2)*pinv(D_T(2,2)), D_T(1,1)-D_T(1,2)*pinv(D_T(2,2))*D_T(2,1); pinv(D_T(2,2)), -pinv(D_T(2,2))*D_T(2,1)];

% Scaling: z-->z/rho, out-->out*K
A_L = rho*A_L;
B_L = sqrt(rho)*B_L;
C_L(1,:) = C_L(1,:)*K;
C_L = C_L*sqrt(rho);
D_L(1,:) = D_L(1,:)*K;

if cnjgt==true
    % Now we drive the conjugate system: z-->1/z
    A_L = pinv(A_L');
    L = ss(A_L, -A_L*C_L', B_L'*A_L, D_L'-B_L'*A_L*C_L', 1);
else
    L = ss(A_L, B_L, C_L, D_L, 1);
end

L = tf(L);
L11 = tf(real(L(1,1).num{1}), real(L(1,1).den{1}), -1);
L12 = tf(real(L(1,2).num{1}), real(L(1,2).den{1}), -1);
L21 = tf(real(L(2,1).num{1}), real(L(2,1).den{1}), -1);
L22 = tf(real(L(2,2).num{1}), real(L(2,2).den{1}), -1);
L = [L11, L12; L21, L22];
L = ss(L);

end