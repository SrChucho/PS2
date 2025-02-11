function [v, sprime, y, yt] = bellman(theta, s, p)

sprime     = golden('valfunc', 1e-10, s - 1e-10, s, theta, p); 

v          = valfunc(sprime, s, theta, p); % evaluates the val fn at optimal choice
    
if nargout > 2

   [v, vt] = valfunc(sprime, s, theta, p);

   y       = funbas(p.fspace, s)*theta -  v;% residual in Bellman equation (want to = 0), Eq (5)

   yt      = funbas(p.fspace, s) - vt; % and its jacobian, Eq (6), derivative wrt theta

end

