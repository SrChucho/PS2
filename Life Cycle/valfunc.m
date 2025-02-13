function v = valfunc(c, a, state, thetae, thetau, p)
    if any(c <= 0)
        v = -Inf;
        return;
    end

    if strcmp(state, 'e')
        income = 1;  % normalized labor income when employed
    else
        income = p.b;  % unemployment benefits
    end
    aprime = (1+p.r)*a + income - c;

    % Check asset accumulation constraint
    if any(aprime < p.amin)  
        v = -Inf;
        return;
    end

    % interpolate value function
    thetaeprime = funbas(p.fspace, aprime)*thetae ;
    
    thetauprime = funbas(p.fspace, aprime)*thetau ; 

    if strcmp(state, 'e')  % employed
        
        v = log(c) + p.beta*((1-p.sigma)*thetaeprime + p.sigma*thetauprime  ) ; % Eq (3)
    
    else                   % unemployed
    
        v = log(c) + p.beta*(p.lambda*thetaeprime + (1 - p.lambda)*thetauprime  )  ; % Eq (5)
    
    end