

function G = Generate_System(ord, K, rho)

warningState = warning('query', 'Control:analysis:NormInfinite3');
warning('off', 'Control:analysis:NormInfinite3');
G = drss(ord); G.D = 0;
while max(abs(pole(G))) > 1/rho || hinfrho(G, rho) > K
    [~, ~] = lastwarn;
    G = drss(ord);
    G.D = 0;
end
warning(warningState);

% z = tf('z', 1);
% G = 0.2*(z-0.1)/((z-0.8)*(z^2+z+0.4));

% % For this system, the reduced model obtained via Hankel singular values
% % trucation becomes unstable.
% G = (-0.1273*z^2 - 0.2194*z - 0.007219) / (z^3 + 1.456*z^2 + 1.067*z + 0.2734);
% % rng(3), N=Nt=7, SNR=10, order=3

% G = (0.0045*z - 0.007916) / (z^3 - 0.3306*z^2 - 0.1487*z - 0.1659);
% rng(3), N=Nt=20, SNR=10, order=3

end