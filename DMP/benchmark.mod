%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
var S lambdaw lambdaf kappa z theta;      % 5 endogenous variables
    
varexo ez; %                        % 1 exogenous variable

parameters beta  sigma  gamma  b  rhoz  sz  alpha  B; % 8 parameters

load dynareinput; % read inputs

% Assigned parameters
beta  = xpar(1); 
sigma = xpar(2);
gamma = xpar(3);
b     = xpar(4);
rhoz  = xpar(5);
sz    = xpar(6);
alpha = xpar(7);
B     = xpar(8);
%kappa = xpar(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODEL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;                              % 5 equations

% total surplus equation, from Nash bargaining
S       = z - b + (1 - sigma - gamma*lambdaw)*beta*S(+1); % Eq (15) in PS2 document

% AR proc
%logz = log(z);
log(z)  = rhoz*log(z(-1)) + sz*ez;                       % Eq (1) in PS2 documen

% worker's matching rate                                  % Eq (18) in PS2 document
lambdaw = B^(1/alpha)*((1-gamma)/kappa)^((1-alpha)/alpha)*(beta*S(+1))^((1-alpha))/alpha; 
%lambdaw  = B ;
//% * theta^(1-alpha);

% firms matching rate, from lambdaf = m/v and theta = v/u
lambdaf = lambdaw/theta;
%lambdaf  = B ;
//% *theta^(-alpha);

% cost of posting a vacancy equal to PV of profits from employing
% Thus, free entry condition is eq (8) kappa = lambdaf * beta * E(J(z'))
% note that from eq (7) J(z) = (1-gamma)*S(z), thus
% by free condition in eq (16) kappa = B*theta^-alpha*(1-gamma)*beta*E(z')
% since theta = 1, 
kappa   = B*(1-gamma)*beta*S(+1);                      

% theta, eq (17) 
theta = ((1 - gamma)* B/kappa)^(1/alpha)*(beta*S(+1))^(1/alpha);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END MODEL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%


shocks;  %

var ez = 1; 

end;

initval;
S       = xsteady(1);
lambdaw = xsteady(2);
lambdaf = xsteady(3);
theta   = xsteady(4);
kappa   = xsteady(5);
z       = 1;
end;

resid; 
steady;

stoch_simul(order = 1, nograph);