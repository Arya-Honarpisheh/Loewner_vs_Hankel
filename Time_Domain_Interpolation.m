close all
clear
clc

%%
seed = 0;
rng(seed)
set(groot, ['Default', 'Stair', 'LineWidth'], 2);

% G is the true system
% Gci is the central intepolant
% Gro is the rank minimized system

K = 1; rho = 1/0.95; % Stability Margins
Nt = 60; % We have Nt time-domain samples
N = 30; % The number of Markov parameters used in the identification proccess
SNR = 20; % signal to noise ratio
R = diag(rho.^(0:N-1));
Ri = diag(rho.^(-(0:N-1)));
R2 = diag(rho.^(2*(0:N-1)));
R2i = diag(rho.^(-2*(0:N-1)));

%%%%%%%%%%% Generating the true system %%%%%%%%%%%%%%%%%%
G = Generate_System(8, K, rho);

%%%%%%%%%%%% Generating the training data %%%%%%%%%%%%%%%
[u, y, Tu, epsilon_t] = Generate_TrainingData(G, Nt, SNR);
%Tu is Nt*Nt

%%

Options = sdpsettings('solver', 'mosek', 'verbose',0, 'debug',1);

gci = sdpvar(N, 1, 'full', 'real'); % gci is the central interpolant
Tgci = toeplitz(gci, [gci(1), zeros(1,N-1)]);
Tu_N = toeplitz(u, [u(1), zeros(1,N-1)]); % Tu_N is Nt*N
assign(gci, Tu_N\y);
%%%%%%%%%%%%% Specifying the Consistency Set %%%%%%%%%%%%%%%%
ConstraintsHankel = [1/K*norm(R*Tgci*Ri, 2) <= 1, norm(y-Tu_N*gci, 'inf') <= epsilon_t];
%%%%%% Finding an impulse response in the consistency set %%%%%%
diagnostics = optimize(ConstraintsHankel,[],Options);
if diagnostics.problem == 1
    error("The data is inconsistent with the apriori information.")
end
% now we can use gci to parametrize all consistent interpolants

gci = value(gci);
Cmin = [1, zeros(1,N-1)];
Cpos = gci'*R / K;
Api = [zeros(N-1,1), eye(N-1); 0, zeros(1,N-1)];

%%%%%%%%%% Finding LFT for the parametrization %%%%%%%%%%%%%%%%
L = CminCposApi2LFT(Cmin, Cpos, Api, rho, K, true);
clc

%%
disp("___________________________________________________________________")
[h, p, v] = HankelRankMinimization(L, N, y, Tu_N, epsilon_t, rho);
%%%%%%%%%%%%%%%%%%%%%%%%% Finding Q %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cmin = v'*R;
Cpos = p'*R;
% Api = Api;
LQ = CminCposApi2LFT(Cmin, Cpos, Api, rho, 1, true);
Q = minreal(LQ(1,1), [], false);
GidH = lft(L, Q); % implementing minreal here may lead to a change in
                  % the Hinfnorm
% This is one of the GidH we can get. Notice that using a different Q
% satisfying (Q(v) = p in BH_inf,1) results in a different GidH.
% Nonetheless, their impulse responses should be closed to each other
%%%%%%%%%%%%%%% Obtaining the reduced model %%%%%%%%%%%%%%%%%%%
if isempty(GidH.A)
    % The system is identified to be identically zero with a specified
    % initial condition
    GroH = GidH;
else
    GroH = Hankel_Reduction(h, 0.9999);
    % There is no gaurantee that the system stays in the BH_inf(K,rho) once
    % we do the balanced reduction, still as we minimized the rank
    % of the system, it should be really closed to a low-rank system.
    % Therefore, there is a great chance that the reduced order system
    % behaves similar to the identified one.
end

%%

[za, wa, h, p, v, dcg] = LoewnerRankMinimization(L, N, y, Tu_N, epsilon_t, rho);
%%%%%%%%%%%%%%%%%%%%%%%%% Finding Q %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Obtaining the reduced model %%%%%%%%%%%%%%%%%%%
    GroL = Loewner_Reduction(wa, conj(wa), za, conj(za), dcg, 0.9999);
    % There is no gaurantee that the system stays in the BH_inf(K,rho) once
    % we do the balanced reduction, still as we minimized the rank
    % of the system, it should be really closed to a low-rank system.
    % Therefore, there is a great chance that the reduced order system
    % behaves similar to the identified one.

%%

n = 30; Nmax = 2*n;
theta = linspace(pi/(n+1), pi-pi/(n+1), n)';
za = exp(theta*1i) / rho;
zb = conj(za);
L = zeros(n,n);
Lro = zeros(n,n);
for i = 1:n
    for j = 1:n
        L(i,j) = (evalfr(G, zb(i)) - evalfr(G, za(j))) / (zb(i) - za(j));
        Lro(i,j) = (evalfr(GroL, zb(i)) - evalfr(GroL, za(j))) / (zb(i) - za(j));
    end
end

[g, ~] = impulse(G, Nmax-1);
[gidH, ~] = impulse(GidH, Nmax-1);
[groH, ~] = impulse(GroH, Nmax-1);
[groL, ~] = impulse(GroL, Nmax-1);

nH = floor(Nmax/2); mH = Nmax - nH;

H = hankel(g(2:nH+1)', g(nH+1:end));
Hid = hankel(gidH(2:nH+1)', gidH(nH+1:end));
Hro = hankel(groH(2:nH+1)', groH(nH+1:end));

%%%%%%%%%%%%%%%%%%%%%%%% results %%%%%%%%%%%%%%%%%%%%%%

set(groot, ['Default', 'Line', 'LineWidth'], 2)
% set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.1)
set(groot, 'DefaultAxesFontSize', 24)
set(groot, 'DefaultAxesFontWeight', 'bold')

figure; hold on
stairs(0:Nmax, [g; g(end)])
stairs(0:Nmax, [groH; groH(end)])
stairs(0:Nmax, [groL; groL(end)])
% stairs(0:Nmax, [gidH; gidH(end)])
legend("True", "HRM", "LRM") % , "IdH")
title("Impulse Responses", FontSize=24)
xlim([0, 40]);

figure; hold on
subplot(1,2,1); hold on
St = svd(H)/sum(svd(H)); Sr = svd(Hro)/sum(svd(Hro));
bar(1:15, [St(1:15), Sr(1:15)])
legend("True", "HRM")
title("Hankel Singular Values", FontSize=24)
subplot(1,2,2); hold on
St = svd(L)/sum(svd(L)); Sr = svd(Lro)/sum(svd(Lro));
bar(1:15, [St(1:15), Sr(1:15)])
legend("True", "LRM")
title("Loewner Singular Values", FontSize=24)

disp("___________________________________________________________________")
show_sys(G, 'True', rho)
show_sys(GidH, 'High-order Identified', rho)
show_sys(GroH, 'Reduced-Order (Hankel)', rho)
show_sys(GroL, 'Reduced-Order (Loewner)', rho)

fprintf("Hinf norm error for identified (Hankel) system is: %f \n", norm(G-GidH, inf))
fprintf("Hinf norm error for reduced-order (Hankel) system is: %f \n", norm(G-GroH, inf))
fprintf("Hinf norm error for reduced-order (Loewner) system is: %f \n", norm(G-GroL, inf))

figure; hold on
%%%%
subplot(2,2,1);
pzmap(G)
title("True System")
axis([-1,1,-1,1])
%%%%%
subplot(2,2,2);
pzmap(GidH)
title("High-order Identified System")
axis([-1,1,-1,1])
%%%%%
subplot(1,2,1);
pzmap(GroH)
title("Identified System using HRM")
axis([-1,1,-1,1])
%%%%%
subplot(1,2,2);
pzmap(GroL)
title("Identified System using LRM")
axis([-1,1,-1,1])

figure; hold on
[mag,~,wout] = bode(G);
Mag=20*log10(mag(:)); wout = log10(wout);
loglog(wout,Mag); grid on;
[mag,~,wout] = bode(GroH);
Mag=20*log10(mag(:)); wout = log10(wout);
loglog(wout,Mag); grid on;
[mag,~,wout] = bode(GroL);
Mag=20*log10(mag(:)); wout = log10(wout);
loglog(wout,Mag); grid on;
% bode(GidH)
legend("True", "HRM", "LRM") % , "IdH")
title("Frequency Responses", FontSize=24)
xlabel("Frequency")
ylabel("Magnitude")



























