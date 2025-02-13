clear;
close all;
clc;

set(groot, 'DefaultAxesLineWidth', 1.5);
set(groot, 'DefaultLineLineWidth', 3);
set(groot, 'DefaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize',22);

p.beta      = 0.95; 
p.sigma     = 0.10; 
p.lambda    = 0.25; 
p.b         = 0.25; % b < 1
p.r         = 0.05; % updated from PS2 doc   

T           = 50; 

% Quality of Approximation

na          = 501;                                              % number of nodes for a

% Grid of points at which to solve the problem

amin        = 0;                                        
amax        = 10;     
p.amin = amin;  

anode       = amin + (amax - amin)*linspace(0, 1, na)'.^2;      % more points in low a region where more curvature


p.fspace    = fundef({'spli', anode, 0, 3});                    % fn space with cubic splines
a           = funnode(p.fspace);                                % nodes from fn space

na          = numel(a);                                         % number of nodes
Phi         = funbas(p.fspace, a);                              % basis fn at nodes

thetae      = zeros(na, T);                                     % value fn of employed
thetau      = zeros(na, T);                                     % value fn of unemployed

thetace     = zeros(na, T);                                     % consumption fn employed
thetacu     = zeros(na, T);                                     % consumption fn unemployed

thetae(:, T)   = log((1 + p.r)*a + 1) ;                         % terminal value fn,  Eq (7), use that it is optimal to consume all wealth at T to find terminal values
thetau(:, T)   = log((1 + p.r)*a + p.b) ;                       % terminal value fn, Eq (8),  

thetace(:, T)  = (1 + p.r)*a + 1;                               % terminal consumption for employed
thetacu(:, T)  = (1 + p.r)*a + p.b;                             % terminal consumption for unemployed

tol = 1e-6;  % Add this before the loop starts

for t = T - 1 : -1 :  1                                         % iterate backward

disp(t)

Ce          = golden('valfunc', 0.01*ones(na, 1), (1 + p.r)*a + 1, a, 'e', ...
                        thetae(:, t+1), thetau(:, t+1), p); 

E           = valfunc(Ce, a, 'e', thetae(:, t+1), thetau(:, t+1), p);

Cu          = golden('valfunc', 0.01*ones(na, 1), (1 + p.r)*a + p.b, a, 'u', ...
                        thetae(:, t+1), thetau(:, t+1), p); 

U           = valfunc(Cu, a, 'u', thetae(:, t+1), thetau(:, t+1), p);

thetae(:, t)  = E ;                   
thetau(:, t)  = U ;

thetace(:, t) = Ce ;
thetacu(:, t) = Cu ;

end

% plot decision rules and check Euler equation at fine grid of a 

t           = 1;

aa          = linspace(0, 10, 1000)';
Phi         = funbas(p.fspace, aa)  ;

ce          =  Phi * thetace(:,t); 
cu          =  Phi * thetacu(:,t);

ceprime     = Phi*thetace(:,t+1); 
cuprime     = Phi*thetacu(:,t+1);

aprimee     = (1 + p.r)*aa + 1   - ce; 
aprimeu     = (1 + p.r)*aa + p.b - cu;

figure(1)

subplot(2, 2, 1)

plot(aa, [ce, cu])
title('$c$', 'Interpreter', 'Latex')
xlabel('$a$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
legend('employed', 'unemployed')

subplot(2, 2, 2)

plot(aa, [aprimee, aprimeu])
title('$a^{\prime}$', 'Interpreter', 'Latex')
xlabel('$a$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')


subplot(2, 2, 3)

plot(aa, 1/ce - p.beta*(1+p.r)*((1-p.sigma)*(1/(ceprime)+p.sigma*(1/cuprime))) )
title('euler error employed', 'Interpreter', 'Latex')
xlabel('$a$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')

subplot(2, 2, 4)

plot(aa, 1/cu - p.beta*(1+p.r)*((1-p.sigma)*(1/(ceprime)+p.sigma*(1/cuprime))) )
%xlim([min(aa(aprimeu > 1e-3)), 10])
title('euler error un employed', 'Interpreter', 'Latex')
xlabel('$a$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')


% simulate an agent

at = zeros(T, 1);
ct = zeros(T, 1);
yt = zeros(T, 1);

et =  ones(T, 1);                          % 1 = employed

et(10 : 15) = 0; 


for t = 1 : T

    ct(t) = ce(t);                         % use consumption function recovered above
    yt(t) = 1*et(t) + p.b*(1 - et(t));

    if t < T

       at(t+1) =  (1 + p.r)*at(t) + yt(t) - ct(t);

    end
end


figure(2)

subplot(1, 2, 1)

plot([yt, ct])
title('consumption and income', 'Interpreter', 'Latex')
xlabel('$t$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
legend('income', 'consumption')

subplot(1, 2, 2)

plot(at)
title('savings', 'Interpreter', 'Latex')
xlabel('$t$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')



% Calculate MPC for employed state
mpc_e = zeros(na, T);
for t = 1:T
    % Calculate derivative using central differences for interior points
    for i = 2:na-1
        dcda = (thetace(i+1,t) - thetace(i-1,t)) / (a(i+1) - a(i-1));
        mpc_e(i,t) = dcda / (1 + p.r);
    end
    
    % Forward difference for first point
    dcda = (thetace(2,t) - thetace(1,t)) / (a(2) - a(1));
    mpc_e(1,t) = dcda / (1 + p.r);
    
    % Backward difference for last point
    dcda = (thetace(na,t) - thetace(na-1,t)) / (a(na) - a(na-1));
    mpc_e(na,t) = dcda / (1 + p.r);
end


figure('Position', [100 100 800 500])
t_plot = 1; % 
plot(a, mpc_e(:,t_plot), 'b-', 'LineWidth', 2)
xlabel('Assets')
ylabel('MPC')
title(['MPC for employed= ' num2str(t_plot)])
grid on


% Calculate MPC for employed state
mpc_e = zeros(na, T);
for t = 1:T
    % Calculate derivative 
    for i = 2:na-1
        dcda = (thetace(i+1,t) - thetace(i-1,t)) / (a(i+1) - a(i-1));
        mpc_e(i,t) = dcda / (1 + p.r);
    end
    
    % Forward difference for first point
    dcda = (thetace(2,t) - thetace(1,t)) / (a(2) - a(1));
    mpc_e(1,t) = dcda / (1 + p.r);
    
    % Backward difference for last point
    dcda = (thetace(na,t) - thetace(na-1,t)) / (a(na) - a(na-1));
    mpc_e(na,t) = dcda / (1 + p.r);
end


t_plot = 1; % 
plot(a, mpc_e(:,t_plot), 'b-', 'LineWidth', 2)
xlabel('Assets')
ylabel('MPC')
title('MPC for employed')
grid on

% Calculate MPC for employed state
mpc_u = zeros(na, T);
for t = 1:T
    % Calculate derivative 
    for i = 2:na-1
        dcda = (thetacu(i+1,t) - thetacu(i-1,t)) / (a(i+1) - a(i-1));
        mpc_u(i,t) = dcda / (1 + p.r);
    end
    
    % Forward difference for first point
    dcda = (thetacu(2,t) - thetacu(1,t)) / (a(2) - a(1));
    mpc_u(1,t) = dcda / (1 + p.r);
    
    % Backward difference for last point
    dcda = (thetacu(na,t) - thetacu(na-1,t)) / (a(na) - a(na-1));
    mpc_u(na,t) = dcda / (1 + p.r);
end


t_plot = 1; % 
plot(a, mpc_u(:,t_plot), 'b-', 'LineWidth', 2)
xlabel('Assets')
ylabel('MPC')
title('MPC for unemployed')
grid on
