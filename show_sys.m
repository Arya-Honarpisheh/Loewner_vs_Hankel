
function [] = show_sys(G, name, rho)
% 
% fprintf("The %s sytem is: \n", name)
% display(tf(G))
fprintf("For %s \n", name)
fprintf("maximum absolute values of the poles: %f \n", max(abs(pole(G))))
G = ss(G);
fprintf("order of the system: %d \n", length(G.A))
fprintf("H_inf,rho norm: %f \n", hinfrho(G, rho))
disp("___________________________________________________________________")

end