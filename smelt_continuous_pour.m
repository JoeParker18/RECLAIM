%% SMELT continuous pour — four-state coupled model (corrected energy split)
% Units: ALL temperatures in degrees Celsius.
clear; clc;
%% ---- Parameters ----
P     = 3000;        % W
eta   = 0.45;        % provisional coupling efficiency
m_Al_total = 0.200;  % kg, total charge
m_c        = 0.358;  % kg, crucible
T0    = 20;  Tinf = 20;  Tm = 660;  Tpour = 700;   % C
Lf     = 397000;     % J/kg
cp_All = 1180;       % J/(kg*C), liquid Al
t_SiC  = 0.01;  k_SiC = 15;  D_i = 0.075;
A_pool = pi*D_i^2/4;
A_int  = 0.0155;
R_int  = t_SiC/(k_SiC*A_int);
UA_loss = 0.337;     % W/C
d_orifice = 0.008;  A_orifice = pi*d_orifice^2/4;
Cd  = 0.80;  g = 9.81;  rho = 2375;
% Way A: charge all solid, at melting point
y0 = [m_Al_total; 0; Tm; Tm];
tspan = [0 3000];
%% ---- Solve ----
odeFun = @(t,y) contPourODE(t, y, m_c, eta, P, R_int, UA_loss, Tinf, Tm, ...
                            Lf, cp_All, Cd, A_orifice, A_pool, g, rho);
opts   = odeset('Events', @depletionEvent, 'MaxStep', 1.0, 'RelTol', 1e-6);
[t, y, te, ye] = ode45(odeFun, tspan, y0, opts);
%% ---- Results ----
disp(['Continuous pour duration = ', num2str(te), ' s  (', num2str(te/60), ' min)']);
%% ---- Plots ----
figure;
subplot(2,1,1);
plot(t, y(:,1), 'c-', t, y(:,2), 'm-', 'LineWidth', 1.8);   % cyan solid, yellow liquid
xlabel('Time (s)', 'FontSize', 16); ylabel('Mass (kg)', 'FontSize', 16);
legend('Solid m_{solid}', 'Liquid pool m_{liquid}', 'FontSize', 14);
set(gca, 'FontSize', 16); grid on;
subplot(2,1,2);
plot(t, y(:,3), 'r-', t, y(:,4), 'g-', 'LineWidth', 1.8);   % red crucible, green pool
xlabel('Time (s)', 'FontSize', 16); ylabel('Temperature (\circC)', 'FontSize', 16);
legend('Crucible T_c', 'Pool T_{liquid}', 'Location','southeast', 'FontSize', 14);
set(gca, 'FontSize', 16); grid on;
%% ---- Functions ----
function dydt = contPourODE(t, y, m_c, eta, P, R_int, UA_loss, Tinf, Tm, ...
                            Lf, cp_All, Cd, A_orifice, A_pool, g, rho)
    m_solid  = max(y(1),0);
    m_liquid = max(y(2),0);
    T_c      = y(3);
    T_liquid = y(4);
    h = max((m_liquid/rho)/A_pool, 0);
if (m_liquid > 1e-9) && (T_liquid >= Tm)
        mdot_drain = rho * Cd * A_orifice * sqrt(2*g*h);
else
        mdot_drain = 0;
end
% --- energy split: melting first (driven by T_c above Tm), leftover superheats pool ---
    Q_total = max((T_c - Tm)/R_int, 0);
if m_solid > 0 && Q_total > 0
        mdot_melt = Q_total / Lf;
else
        mdot_melt = 0;
end
    Q_leftover = Q_total - mdot_melt*Lf;
    dm_solid  = -mdot_melt;
    dm_liquid = mdot_melt - mdot_drain;
    cp_SiC = 697.088 + 0.656*T_c;   % intercept shifted from 518 (518 + 0.656*273), T_c in C
    Q_loss = UA_loss * (T_c - Tinf);
    dT_c   = (eta*P - Q_total - Q_loss) / (m_c * cp_SiC);
    m_liquid_min = 1e-2;
if m_liquid > m_liquid_min
        dT_liquid = (Q_leftover + mdot_melt*cp_All*(Tm - T_liquid)) / (m_liquid*cp_All);
else
        dT_liquid = 0;
end
    dydt = [dm_solid; dm_liquid; dT_c; dT_liquid];
end
function [value, isterminal, direction] = depletionEvent(t, y)
    value      = (y(1)+y(2)) - 1e-4;
    isterminal = 1;
    direction  = -1;
end