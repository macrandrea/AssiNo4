maxit=1000;
tol=1e-8;
mat=zeros(maxit+1,2);
mat(1,:)=[0,0];
for i=1:maxit
    x=mat(i,1);
    y=mat(i,2);
    H=[12*x^2, 1;
       1 , 2];
    J=[4*x^3+y; x+2*y+2];
    p=H\J;
    mat(i+1,:) = mat(i,:)-p';
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end   
end