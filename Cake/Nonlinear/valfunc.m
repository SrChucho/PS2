function [v, vt] = valfunc(sprime, s, theta, p)

c = s - sprime;

Phinode = funbas(p.fspace,sprime);

v          = exp((1-p.beta)*log(c) + p.beta*log(Phinode*theta)) ;   % RHS of Bellman equation (4)

 % derivative w.r.t theta coefficients (1st coefficient is in 1s column, 2nd in 2nd column ...)
vt         = p.beta*exp((1-p.beta)*log(c) + p.beta*log(Phinode*theta)).*(1./(Phinode*theta)).*Phinode   ;  

