function [y, yx] = solvex(x, p)

Si       = x(1        : p.nz);
lambdawi = x(p.nz + 1 :  end);

y     = zeros(p.nz*2, 1);

y(1        : p.nz)  = Si - (p.zi - p.b + p.beta*(1 - p.sigma - p.gamma*lambdawi).*p.Pi*Si); % this is  Eq (15)
y(p.nz + 1 :  end)  =  lambdawi - p.B^(1/p.alpha)*(((1-p.gamma)/p.kappa)*(p.beta*p.Pi*Si)).^((1-p.alpha)/p.alpha); % this is  Eq (18)
    
yx       = zeros(p.nz*2, p.nz*2); % jacobian

yx(1 : p.nz,        1 :   p.nz) = eye(p.nz) - p.beta*(1 - p.sigma - p.gamma*lambdawi).*p.Pi; % derivative of 1st eq wrt to Si
yx(1 : p.nz, p.nz + 1 : 2*p.nz) =  p.beta*(p.gamma).*p.Pi*Si.*eye(p.nz) ; % derivative of 2nd eq wrt to lambdawi

yx(p.nz + 1 : 2*p.nz,        1 :   p.nz) = -p.B^(1/p.alpha)*((1-p.gamma)/p.kappa)^((1-p.alpha)/p.alpha)*((1-p.alpha)/p.alpha)*(p.beta*p.Pi).^((1-p.alpha)/p.alpha)*Si.^((1-2*p.alpha)/p.alpha)*ones(1, p.nz); %derivative of second equation wrt to Si 
yx(p.nz + 1 : 2*p.nz, p.nz + 1 : 2*p.nz) = eye(p.nz);; % derivative of second equation wrt to lambdawi

end