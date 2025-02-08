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
p.b         = 0.25;
p.r         = 0.04;   

T           = 50; 

% Quality of Approximation

na          = 501;                                              % number of nodes for a

% Grid of points at which to solve the problem

amin        = 0;                                        
amax        = 10;     

anode       = amin + (amax - amin)*linspace(0, 1, na)'.^2;      % more points in low a region where more curvature

p.fspace    = fundef({'spli', anode, 0, 3});
a           = funnode(p.fspace);

na          = numel(a);
Phi         = funbas(p.fspace, a);

thetae      = zeros(na, T);                                     % value function
thetau      = zeros(na, T); 

thetace     = zeros(na, T);                                     % consumption function
thetacu     = zeros(na, T); 

thetae(:, T)   = ...;                                           % use that it is optimal to consume all wealth at T to find terminal values
thetau(:, T)   = ...;

thetace(:, T)  = ...;
thetacu(:, T)  = ...;


for t = T - 1 : -1 :  1                                         % iterate backward

disp(t)

Ce          = golden('valfunc', 0.01*ones(na, 1), (1 + p.r)*a + 1, a, 'e', thetae(:, t+1), thetau(:, t+1), p); 

E           = valfunc(Ce, a, 'e', thetae(:, t+1), thetau(:, t+1), p);

Cu          = golden('valfunc', 0.01*ones(na, 1), (1 + p.r)*a + p.b, a, 'u', thetae(:, t+1), thetau(:, t+1), p); 

U           = valfunc(Cu, a, 'u', thetae(:, t+1), thetau(:, t+1), p);

thetae(:, t)  = ...;                   
thetau(:, t)  = ...;

thetace(:, t) = ...;
thetacu(:, t) = ...;

end


% plot decision rules and check Euler equation at fine grid of a 

t           = 1;

aa          = linspace(0, 10, 1000)';

ce          = ....; 
cu          = ...;

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

plot(aa, ...)
title('euler error employed', 'Interpreter', 'Latex')
xlabel('$a$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')

subplot(2, 2, 4)

plot(aa, ...)
xlim([min(aa(aprimeu > 1e-3)), 10])
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

    ct(t) = ...;                           % use consumption function recovered above
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

