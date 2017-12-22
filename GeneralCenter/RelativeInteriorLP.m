function [xi,x0,lambda_I] = RelativeInteriorLP(A,b,isEquality,isInequality)

[~,d] = size(A);  

f = -[ zeros(1,d) 1];

A_eq = [A(isEquality,:) zeros(sum(isEquality),1)];
b_eq = b(isEquality);

A_ineq = [A(isInequality,:) ones(sum(isInequality),1)];
b_ineq = b(isInequality);

upper_bound = ones(d+1,1)*inf;
upper_bound(end) = 1;

options = optimoptions('linprog','Algorithm','dual-simplex','Display','none');
[x0,~,exitflag,output,lambda] = linprog(f,A_ineq,b_ineq,A_eq,b_eq,[],upper_bound,options);
assert(exitflag == 1,'Error using linprog in ChebyshevCenter');

lambda_I = lambda.ineqlin;

xi = x0(end);
x0 = x0(1:end-1);

end