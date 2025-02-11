clear;
close all;
clc;

set(groot, 'DefaultAxesLineWidth', 1.5);
set(groot, 'DefaultLineLineWidth', 3);
set(groot, 'DefaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize',22);


p.beta      = 0.95; 

% Quality of Approximation

smin        = 1e-9; 
smax        = 5; 

ns          = 251;                    % number of nodes for s

switch 'cheb'

    case 'cheb'

        p.fspace    = fundefn('cheb', ns, smin, smax);
        s           = funnode(p.fspace);

    case 'spli'

        snode       = nodeunif(ns - 2, 0, 1).^3*(smax - smin) + smin; % cubic transformation
        p.fspace    = fundef({'spli', snode, 0, 3}); % defines function APROX space 
        s           = funnode(p.fspace);             % nodes

end

Phi         = funbas(p.fspace, s);   % evaluate matrix of nxn functions at n points in s

v           = zeros(ns, 1);          % guess initial value
theta       = Phi\v;                 % unknown coefficients


for i = 1 : 100

    thetaold   = theta;

    sprime     = golden('valfunc', 1e-10, s - 1e-10, s, theta, p);        % edit valfunc
    
    c = s - sprime;

    theta      = (funbas(p.fspace, s) - p.beta*funbas(p.fspace, sprime))\log(c) ; % Eq (8)

    err        = norm(theta - thetaold)/norm(theta);

    fprintf('%4i %6.2e \n', [i, err]);    


    if err < 1e-9 , break, end

end

s          = nodeunif(10000, smin, smax);
sprime     = golden('valfunc', 1e-10, s - 1e-10, s, theta, p); 

v          =  funbas(p.fspace, s)*theta;   % LHS of Bellman 

sprimetrue = p.beta*s;  % This is Eq (5)

vtrue      = log(1 - p.beta)/(1 - p.beta) + p.beta/(1 - p.beta)^2*log(p.beta) + 1/(1 - p.beta)*log(s);



% plot decision rules and value

figure(1)

subplot(2, 2, 1)

plot(s, sprime)
hold on
plot(s, sprimetrue, '--')
title('$s^{\prime}$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
legend('numerical', 'actual')
xlim([0.01, 5])

subplot(2, 2, 2)

plot(s, v)
hold on
plot(s, vtrue, '--')
title('$v$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
xlim([0.01, 5])


subplot(2, 2, 3)

plot(s, sprime - sprimetrue)

title('error $s^{\prime}$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
xlim([0.01, 5])


subplot(2, 2, 4)

plot(s, v - vtrue)
title('error $v$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
xlim([0.01, 5])

% residuals

figure(2)

plot(s, v - funeval(theta, p.fspace, s))
%plot(s, v - valfunc(sprime, s, theta, p))

title('residuals', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
