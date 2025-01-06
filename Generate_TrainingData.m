

function [u, y, Tu, epsilon_t] = Generate_TrainingData(G, Nt, SNR)

[g, ~] = impulse(G, Nt-1);

u = idinput(Nt, 'rgs') * 1;
Tu = toeplitz(u, [u(1), zeros(1,Nt-1)]);
% yt = lsim(G,u);
yt = Tu*g;
epsilon_t = range(yt)/SNR;
noise = idinput(Nt, 'rbs') * epsilon_t;
y = yt + noise;

figure; hold on
stairs(0:Nt, [yt; yt(end)])
stairs(0:Nt, [y; y(end)])
legend("True", "Noisy")
title("The time-domain data")

end