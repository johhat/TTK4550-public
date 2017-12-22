function [isEquality,isInequality] = AddLinearlyDependent(A,b,isEquality,isInequality)
m_E = sum(isEquality > 0);
[m,~] = size(A);

if m_E == 0 || m_E == m
    return
end

A_E = A(isEquality,:); 
b_E = b(isEquality);
A_I = A(isInequality,:); 
b_I = b(isInequality);

A_aug = [A_E' A_I' ; b_E' b_I'];
A_rrechelon = rref(A_aug);

r = rank(A_rrechelon(:,1:m_E));

for i = m_E+1:m
    % if column i is a linear combination of the first m_E columns, add to I
    c = A_rrechelon(r+1:end,i);
    if all(c == 0)
        isEquality(i) = true;
        isInequality(i) = false;
    end
end

end