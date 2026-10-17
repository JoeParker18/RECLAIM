%% SMELT two-node batch melt — parameters & derived quantities
% Units: ALL temperatures in degrees Celsius.
clear; clc;

%% ---- Power / efficiency ----
P     = 1000;        % W, input power
eta   = 0.45;        % provisional coupling efficiency (report benchmark)

%% ---- Masses ----
m_Al  = 0.200;       % kg, aluminum charge
m_c   = 0.358;       % kg, measured crucible mass

%% ---- Temperatures (C) ----
T0    = 20;          % start / ambient temp
Tinf  = T0;          % surroundings temp for loss term
Tm    = 660;         % Al melting point (660 C)
Tpour = 700;         % target pour temp (stop condition), ~710-810 C band

%% ---- Latent heat ----
Lf    = 397000;      % J/kg, Al fusion
dTpc  = 5;           % phase-change smoothing window, ~5-10 C (interval -> unchanged)
% Geometry
d_orifice = 0.008;              % m, orifice diameter (8 mm)
D_i       = 0.075;              % m, crucible bore (pool diameter)
A_orifice = pi*d_orifice^2/4;   % m^2
A_pool    = pi*D_i^2/4;         % m^2

% Flow
Cd  = 0.80;                     % discharge coeff, short-tube L/D~5.2
g   = 9.81;                     % m/s^2 (Earth)
rho = 2375;                     % kg/m^3, molten Al
dP  = 0;                        % Pa, purge-and-vent

% Initial condition
h0  = 0.047;                    % m, starting pool height
%% ---- Specific heat fits, J/(kg*C) — VERIFY before citing ----
% cp_Als(T) = 887.12  + 0.44*T    (solid Al,  T in C)
% cp_All    = 1180                (liquid Al, ~constant)
% cp_SiC(T) = 697.088 + 0.656*T   (SiC,       T in C)
% (intercepts shifted from the Kelvin fits so cp is identical at the same
%  physical temperature: 767 + 0.44*273 = 887.12; 518 + 0.656*273 = 697.088)
cp_All = 1180;

%% ---- Internal resistance R_int (crucible -> Al) ----
% BN coat dropped as thermally negligible vs SiC wall (~0.4% of SiC term).
% Keep k_SiC LOW: measured 358 g crucible implies high porosity.
t_SiC  = 0.01;                           % m, your 1 cm wall
k_SiC  = 15;                             % W/(m*C), low end for porous cast SiC
D_i    = 0.075;                          % m, crucible inner diameter (7.5 cm)
A_floor = pi*D_i^2/4;                    % m^2, crucible floor area
V_Al   = 0.000208;                       % m^3, molten Al volume (~500 g)
H_fill = V_Al / A_floor;                 % m, fill height
A_int  = pi*D_i*H_fill + A_floor;        % m^2, wetted wall + floor
R_int  = t_SiC/(k_SiC*A_int);
h_stop = 0.001;   % m, residual heel — stop drain here, not at exactly 0

%% ---- Loss conductance (U_loss * A_loss) ----
U_loss  = 2.674;     % W/(m^2*C), your verified hand calc
A_loss  = 0.12585;   % m^2, insulation casing outer area
UA_loss = U_loss * A_loss;   % W/C, conductance -> use THIS in loss term

%% ---- Sanity checks ----
disp(['R_int   = ', num2str(R_int), ' C/W']);
disp(['UA_loss = ', num2str(UA_loss), ' W/C   (expect ~0.337)']);

%% ODE45 solver portion
function dTdt = smeltODE(t, T, eta, P, m_c, m_Al, R_int, UA_loss, Tinf, Tm, dTpc, Lf, cp_All)
    Tc  = T(1);
    Tal = T(2);
    cp_SiC = 697.088 + 0.656*Tc;
    cp_Als = 887.12  + 0.44*Tal;
    if Tal < Tm - dTpc/2
        cp_Al_eff = cp_Als;
    elseif Tal <= Tm + dTpc/2
        cp_Al_eff = cp_Als + Lf/dTpc;
    else
        cp_Al_eff = cp_All;
    end
    Q_cAl  = (Tc - Tal)/R_int;
    Q_loss = UA_loss*(Tc - Tinf);
    dTc  = (eta*P - Q_cAl - Q_loss)/(m_c*cp_SiC);
    dTal = Q_cAl/(m_Al*cp_Al_eff);
    dTdt = [dTc; dTal];
end

smeltODE(0, [T0; T0], eta, P, m_c, m_Al, R_int, UA_loss, Tinf, Tm, dTpc, Lf, cp_All)

%% Step 2 of the ode creation
function [value, isterminal, direction] = pourEvent(t, T, Tpour)
    value = T(2) - Tpour;
    isterminal = 1;
    direction = 1;

end

%% ---- Solve ----
% Wrap the functions so ode45 sees only (t,T); the rest are filled from workspace
odeFun   = @(t,T) smeltODE(t, T, eta, P, m_c, m_Al, R_int, UA_loss, Tinf, Tm, dTpc, Lf, cp_All);
eventFun = @(t,T) pourEvent(t, T, Tpour);

opts   = odeset('Events', eventFun);        % register the stop condition
tspan  = [0 3600];                          % generous ceiling (s); event halts early
T_init = [T0; T0];                          % both nodes start at ambient

[t, T, te, ye] = ode45(odeFun, tspan, T_init, opts);

%% ---- Results ----
melt_time = te;                             % event time = melt duration (s)
disp(['Melt time = ', num2str(melt_time), ' s  (', num2str(melt_time/60), ' min)']);

%% ---- Plot ----
figure;
plot(t, T(:,1), 'r-', t, T(:,2), 'b-', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 16);
ylabel('Temperature (\circC)', 'FontSize', 16);
legend('Crucible T_c', 'Aluminum T_{al}', 'Location', 'southeast', 'FontSize', 16);
set(gca, 'FontSize', 18);

grid on;

%% Solving the time it takes to have the aluminum flow through the hole

function dhdt = drainODE(t, h, Cd, A_orifice, A_pool, g, rho, dP)
    % guard against tiny negative h from solver overshoot
    h_eff = max(h, 0);
    v = Cd * sqrt(2*(g*h_eff + dP/rho));   % orifice exit velocity
    dhdt = -(A_orifice/A_pool) * v;         % pool height drops as metal leaves
end

drainODE(0, h0, Cd, A_orifice, A_pool, g, rho, dP)
eventFun = @(t,h) drainEvent(t, h, h_stop);
odeFun = @(t,h) drainODE(t, h, Cd, A_orifice, A_pool, g, rho, dP);
opts = odeset('Events', eventFun);
[t, h, te, he] = ode45(odeFun, [0 60], h0, opts);

function [value, isterminal, direction] = drainEvent(t, h, h_stop)
    value      = h - h_stop;   % zero when pool drops to the stop height
    isterminal = 1;            % stop integration there
    direction  = -1;           % trigger on the DOWNWARD crossing (h falling)
end
disp(['Drain time = ', num2str(te), ' s']);