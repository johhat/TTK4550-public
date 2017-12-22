function [A_C,b_C,lambda,x_c,r_opt] = ChebyshevCenterLP(A,b)
    [m,d] = size(A);
    f = -[zeros(1,d) 1]; % Objective function
    
    % Normalize constraints by length
    n = vecnorm(A,2,2);
    
    if any(n==0)
       assert(false,'constraints with vector norm zero discovered'); 
    end
    
    A_norm = A./n;
    b_norm = b./n;
    
    % Construct inequalities used in LP
    A_ineq = [A_norm ones(m,1) ; zeros(1,d) -1];
    b_ineq = [b_norm ; 0];
    
    % Optimize
    options = optimoptions('linprog','Algorithm','dual-simplex','Display','none');
    [x,~,exitflag,~,lambda] = linprog(f,A_ineq,b_ineq,[],[],[],[],options);
    assert(exitflag == 1,'Error using linprog in ChebyshevCenter');
    
    x_c = x(1:end-1);
    r_opt = x(end);
    
    A_C = A_norm;
    b_C = b_norm - r_opt*ones(m,1);
    lambda = lambda.ineqlin(1:end-1);
end

function n = vecnorm(A,p,dim)
    assert(p == 2,'Only p=2 supported');
    assert(dim == 2,'Only dim=2 supported');
    n = sqrt(sum(A.^2,2));
end