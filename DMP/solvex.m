function [y, yx] = solvex(x, p)


Si       = x(1        : p.nz);
lambdawi = x(p.nz + 1 :  end);



y     = zeros(p.nz*2, 1);

y(1        : p.nz)  = Si - (p.zi - p.b + p.beta*(1 - p.sigma - p.gamma*lambdawi).*p.Pi*Si); 
y(p.nz + 1 :  end)  = ...;
    

yx       = zeros(p.nz*2, p.nz*2);

yx(1 : p.nz,        1 :   p.nz) = eye(p.nz) - p.beta*(1 - p.sigma - p.gamma*lambdawi).*p.Pi;
yx(1 : p.nz, p.nz + 1 : 2*p.nz) = ...;

yx(p.nz + 1 : 2*p.nz,        1 :   p.nz) = ...;
yx(p.nz + 1 : 2*p.nz, p.nz + 1 : 2*p.nz) = ...;

end