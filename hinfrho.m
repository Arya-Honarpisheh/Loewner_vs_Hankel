
function [hinf ,freq] = hinfrho(G, rho)
% This function computes the Hinf norm of the transfer function G
% on the circle {|z| = 1/rho} assuming that G is holomorphic on
% {|z|>= 1/rho}
% hinf = MAX(abs(G(1/rho*exp(1i*w)))) s.t w \in [0,2pi]
% rho > 1

G = ss(G);
A = G.A * rho; B = G.B * sqrt(rho); C = G.C * sqrt(rho); D = G.D;
Grho = ss(A, B, C, D, 1);
[hinf, freq] = hinfnorm(Grho);

end