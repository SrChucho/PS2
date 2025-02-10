function v = valfunc(sprime, s, theta, p)
Phinode = funbas(p.fspace,sprime);
v  = log(s - sprime) + p.beta*Phinode*theta; % RHS of bellman equation (7) given savings sprime and initial state s and the interpolation coefficients theta
