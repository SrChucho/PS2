function x = golden(f, a, b, varargin)

tol    = 1e-12;

alpha1 = (3-sqrt(5))/2; % 0.3820
alpha2 = (sqrt(5)-1)/2; % 0.6180

d  = b - a; 

x1 = a + alpha1*d; 
x2 = a + alpha2*d; 

f1 = feval(f, x1, varargin{:}); 
f2 = feval(f, x2, varargin{:}); 

d = alpha1*alpha2*d;

x1new = x1;
x2new = x2;
f1new = f1;
f2new = f2;

while any((d) > tol)

 f1 = f1new;
 f2 = f2new;
 x1 = x1new;
 x2 = x2new;
    
 d     = d*alpha2;
 x2new = x1.*(f2<f1) + (x2+d).*(f2>=f1);
 f2new = f1.*(f2<f1) + feval(f, x2 + d, varargin{:}).*(f2>=f1);
  
 x1new = (x1-d).*(f2<f1)+x2.*(f2>=f1);
 f1new = feval(f, x1 - d, varargin{:}).*(f2<f1)+f2.*(f2>=f1);
  
end

x = x2new.*(f2new>=f1new) + x1new.*(f2new<f1new);