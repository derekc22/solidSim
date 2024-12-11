function [tList, P0List] = solidSim(rho, T0, k, n, gamma, Mbar, ri, rf, L, Astar)

Patm = 1.01E5;

init = [Patm, ri];
tspan = [0, 20];
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9, 'Events', @stopAfterBurnout);




[t, f, te, ye, ie] = ode45(@odefun, tspan, init, options);



function [value,isTerminal,direction] = stopAfterBurnout(t, f)
value(1) = f(1) - Patm;
isTerminal(1) = 1;
direction(1) = -1; 
end



%% Plot 

tList = t;
P0List = f(:, 1)/Patm;
plot(tList, P0List, LineWidth=2)

grid on
title("Stagnation Pressure With Respect to Time")
xlabel("Time [s]")
ylabel("Stagnation Pressure [atm]")


%% ODE function

% Let P0 = f(1)
% Let r = f(2)

function dfdt = odefun(t, f)

if f(2) >= rf
    k = 0;
end 


R = 8314/Mbar;

rDot = k*(f(1)/Patm)^n;
V = pi*f(2)^2*L;
mgDot = rDot * 2*pi*f(2) * L * rho;
mnDot = Astar * ((f(1)*sqrt(gamma))/(sqrt(R*T0))) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1))); 

dfdt = [((R*T0)/V)*(mgDot - mnDot);
rDot];


end 
end

