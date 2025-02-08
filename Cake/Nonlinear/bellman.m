function [v, sprime, y, yt] = bellman(theta, s, p)

sprime     = golden('valfunc', 1e-10, s - 1e-10, s, theta, p); 

v          = valfunc(sprime, s, theta, p);
    
if nargout > 2

   [v, vt] = valfunc(sprime, s, theta, p);

   y       = ...;                      % residual in Bellman equation (want to = 0)

   yt      = ...;                      % and its jacobian

end

