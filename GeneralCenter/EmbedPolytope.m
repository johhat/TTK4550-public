function [A_hat,b_hat,V] = EmbedPolytope(A,b,x0)

b_0 = A*x0;
isEquality = (b_0 == b);
isInequality = ~isEquality;

% Basis for nullspace of equality constraints
V = null(A(isEquality,:)); 

A_hat = A(isInequality,:)*V;
b_hat = b(isInequality) - A(isInequality,:)*x0; 

% Remove redundant inequality constraints perpendicular to V
n = sum(A_hat.^2,2);
A_hat = A_hat(n>0,:);
b_hat = b_hat(n>0,:);

end