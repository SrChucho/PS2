clear;
close all;
clc;

set(groot, 'DefaultAxesLineWidth', 1.5);
set(groot, 'DefaultLineLineWidth', 3);
set(groot, 'DefaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize', 22);

darkblue    = [0.02  0.28, 0.46]; % dark blue
maroon      = [0.64, 0.08, 0.18]; % maroon

% Parameters

beta    = 0.96^(1/12);            % discount factor

sigma   = 0.04;                   % probability of separation
gamma   = 0.50;                   % bargaining share worker

b       = 0.7; 
rhoz    = 0.90^(1/12);
sz      = 0.005; 
alpha   = 0.50;                   % elasticity matching function; 

B       = 0.40;    

theta   = 1; 
lambdaf = B;         
lambdaw = B;                      % probability worker matches

% Accuracy

nz      = 51;

zss     = 1;                      % at SS z is 1 

S       = (zss-b)/(1-(1-sigma-gamma*lambdaw)*beta);                    % steady state surplus (20)
kappa   = (1-gamma)*B*(zss-b)/(1/beta - 1 + sigma + gamma*lambdaw);    % fixed cost of vacancy posting (21)

% save steady-state

xsteady = [S; lambdaw; lambdaf; theta];

xpar    = [beta; sigma; gamma; b; rhoz; sz; alpha; B; kappa];

save dynareinput xsteady xpar; 
dynare benchmark noclearall;        % edit this file

% Cleanup

filename = 'benchmark';

eval(['delete ' filename '*.m']);
eval(['delete ' filename '*.mat']);
eval(['delete ' filename '*.swp']);
eval(['delete ' filename '*.log']);
eval(['delete ' filename '*.eps']);
eval(['delete ' filename '*.pdf']);
eval(['delete ' filename '*.fig']);


[zi, Pi]   = rouwenhorst(rhoz, sz, nz);

zi        = exp(zi)';

fprintf('\n')
fprintf('Iterative Method\n')
fprintf('\n')
fprintf('\n')

% Non-linear method 1

lambdawi  = lambdaw*ones(nz, 1);
Si = zeros(nz, 1);
Sj = zeros(nz, 1);

% init 
for iz = 1:nz
    Si(iz) = (zi(iz)-b)/(1-(1-sigma-gamma*B)*beta);
    lambdawi(iz) = B;
end

for i = 1 : 200                                            

    lambdawiold  = lambdawi; 
    ESj = Pi*Si;

    for iz = 1:nz
        Si(iz)       = zi(iz) - b + (1 - sigma - gamma * lambdawi(iz))*beta*ESj(iz) ;         % Eq (22)

        lambdawi(iz) = B^(1/alpha)*((1 - gamma)/kappa)^((1-alpha)/alpha)*(beta*ESj(iz))^((1-alpha)/alpha); % Eq (23)
    end 

    fprintf('%4i %6.2e  \n', [i, norm(lambdawi - lambdawiold)]);    

    if norm(lambdawi - lambdawiold) < 1e-11 , break, end

end



figure(1)
subplot(1, 2, 1)
plot(zi, Si, 'Color', darkblue); % plots data 
hold on
% the next line plots mode prediction
plot(zi, oo_.dr.ys(1) + oo_.dr.ghx(oo_.dr.inv_order_var(1))*(zi - 1), '-.', 'Color', maroon)
set(gca, 'ygrid', 'on')
title('$S(z)$', 'Interpreter', 'Latex')
xlabel('$z$', 'Interpreter', 'Latex')

subplot(1, 2, 2)
plot(zi, lambdawi, 'Color', darkblue); % plots data 
hold on
% the next line plots mode prediction
plot(zi, oo_.dr.ys(2) + oo_.dr.ghx(oo_.dr.inv_order_var(2))*(zi - 1), '-.', 'Color', maroon)
set(gca, 'ygrid', 'on')
title('$\lambda_w(z)$', 'Interpreter', 'Latex')
xlabel('$z$', 'Interpreter', 'Latex')

% Non-Linear Method 2

x0   = [S*ones(nz, 1); lambdaw*ones(nz, 1)]; 

x          = x0; 

fprintf('\n')
fprintf('Newton Method \n')
fprintf('\n')
fprintf('\n')

p.sigma = sigma; 
p.alpha = alpha; 
p.beta  = beta; 
p.B     = B; 
p.gamma = gamma; 
p.kappa = kappa; 
p.b     = b; 
p.zi    = zi; 
p.Pi    = Pi; 
p.nz    = nz; 

for i = 1 : 200                                            

    xold         = x; 

    [y, yx]      = solvex(xold, p);                 % edit this file

    x            = xold - yx\y ;                             % apply Newton's step, x^k+1 = x^k - J(x^k)-1f(x^k)
                                                             % p 167 Judd (1998)

    fprintf('%4i %6.2e %6.2e \n', [i, norm(y), norm(x - xold)]);    

    if norm(y) < 1e-11 && norm(x - xold) < 1e-11 , break, end

end

Sii = x(1        : p.nz);
lambdawii = x(p.nz + 1 :  end);
purple = [0.5, 0, 0.5]; 


figure(1)
subplot(1, 2, 1)
plot(zi, Sii, 'Color', purple, 'Marker', '+'); % plots data 
hold on
% the next line plots mode prediction
plot(zi, oo_.dr.ys(1) + oo_.dr.ghx(oo_.dr.inv_order_var(1))*(zi - 1), '-.', 'Color', maroon)
set(gca, 'ygrid', 'on')
title('$S(z)$', 'Interpreter', 'Latex')
xlabel('$z$', 'Interpreter', 'Latex')

subplot(1, 2, 2)
plot(zi, lambdawii, 'Color', purple, 'Marker', '+'); % plots data 
hold on
% the next line plots mode prediction
plot(zi, oo_.dr.ys(2) + oo_.dr.ghx(oo_.dr.inv_order_var(2))*(zi - 1), '-.', 'Color', maroon)
set(gca, 'ygrid', 'on')
title('$\lambda_w(z)$', 'Interpreter', 'Latex')
xlabel('$z$', 'Interpreter', 'Latex')
