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
p.ns        = ns;
switch 'cheb'

    case 'cheb'

        p.fspace    = fundefn('cheb', ns, smin, smax);
        s           = funnode(p.fspace);

    case 'spli'

        snode       = nodeunif(ns - 2, 0, 1).^3*(smax - smin) + smin; 
        p.fspace    = fundef({'spli', snode, 0, 3});
        s           = funnode(p.fspace);

end

p.Phi       = funbas(p.fspace, s);  % save it to avoid recomputing it 25x25

v           = ones(ns, 1);          % guess initial value
theta       = p.Phi\v;                  % guess theta by projecting v on p.Phi

% a few value function iterations

for i = 1 : 5

    thetaold   = theta;

    [~, ~, y, yt]  = bellman(theta, s, p);   % edit this yourselves

    theta      = thetaold - yt\y; %                     % update theta by projecting v on p.Phi

    err        = norm(theta - thetaold)/norm(theta);

    fprintf('%4i %6.2e \n', [i, err]);    

    if err < 1e-9 , break, end

end


for i = 1 : 100

    thetaold   = theta;

    [~, ~, y, yt]  = bellman(theta, s, p);   % y and yt are the residuals in the Bellman equation and its Jacobian

    theta          = thetaold - yt\y ;  %use y and yt to update theta

    err            = norm(theta - thetaold)/norm(theta);

    fprintf('%4i %6.2e \n', [i, err]);    

    if err < 1e-9 , break, end

end

% Wider grid of nodes to check accuracy

s           = nodeunif(1000, smin, smax);
[v, sprime] = bellman(theta, s, p);


% plot decision rules and value

figure(1)

subplot(2, 2, 1)

plot(s, sprime)
hold on
plot(s, p.beta*s, '--')
title('$s^{\prime}$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')
legend('numerical', 'actual')


subplot(2, 2, 2)

plot(s, v)
hold on
plot(s, (1 - p.beta)*p.beta^(p.beta/(1 - p.beta))*s, '--')
title('$v$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')


subplot(2, 2, 3)

plot(s, sprime - p.beta*s)

title('error $s^{\prime}$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')


subplot(2, 2, 4)

plot(s, v - (1 - p.beta)*p.beta^(p.beta/(1 - p.beta))*s)
title('error $v$', 'Interpreter', 'Latex')
xlabel('$s$', 'Interpreter', 'Latex')
set(gca, 'ygrid', 'on')


% residuals

figure(2)

plot(s, v - funeval(theta, p.fspace, s))

title('residuals', 'Interpreter', 'Latex')

xlabel('$s$', 'Interpreter', 'Latex')

set(gca, 'ygrid', 'on')
