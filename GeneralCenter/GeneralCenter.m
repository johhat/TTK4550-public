function [x_g,r_opt] = GeneralCenter(A,b)
[~,d] = size(A);    
[A_C,b_C,lambda,x_c,r_opt] = ChebyshevCenterLP(A,b);
H = eye(d+1);
[d_C,x0] = PolytopeDimension(A_C,b_C,lambda,x_c);

while d_C > 0
    [A_C,b_C,V] = EmbedPolytope(A_C,b_C,x0);
    H = H*[V x0 ; zeros(size(V,2),1) 1];
    [A_C,b_C,lambda,x_c] = ChebyshevCenterLP(A_C,b_C);
    [d_C,x0] = PolytopeDimension(A_C,b_C,lambda,x_c);
end

x_hom = H*[x0 ; 1];
x_g = x_hom(1:end-1);

end