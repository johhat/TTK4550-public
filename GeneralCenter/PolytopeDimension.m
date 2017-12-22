function [d_P,x0] = PolytopeDimension(A,b,lambda,x0)
    [m,d] = size(A);
    
    % rank(A) <= min(m_E,d). d_N is the rank of the nullspace of A.
    d_N_func = @(A,isEquality) d - rank(A(isEquality,:));
    
    % Arrays of binary values. isI(i) = true iff i is in I, and similarly
    % for E.
    isEquality = lambda > 0;
    isInequality = lambda <= 0;
    
    [isEquality, isInequality] = AddLinearlyDependent(A,b,isEquality,isInequality);
    
    d_N = d_N_func(A,isEquality);
    
    if d_N == 0
        d_P = 0;
        return
    end
    
    while true
        [xi,x0,lambda_I] = RelativeInteriorLP(A,b,isEquality,isInequality);
        
        if xi>0
            d_P = d_N;
            return
        end
        
        lambda = zeros(m,1);
        lambda(isEquality) = 1; % Dummy value
        lambda(isInequality) = lambda_I;
         
        isEquality = lambda > 0;
        isInequality = lambda <= 0;
                 
        [isEquality, isInequality] = AddLinearlyDependent(A,b,isEquality,isInequality);         
        
        d_N = d_N_func(A,isEquality);
         
        if d_N == 0
            d_P = 0;
            return
        end
    end
end