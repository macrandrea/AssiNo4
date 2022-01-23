maxit=10000;
tol=1e-8;
mat=zeros(maxit+1,2);
mat(1,2)=1;
mat(1,1)=1;
for i=1:maxit
    J=[4*(mat(i,1) - 2)^3 + mat(i,2)^2*(2*mat(i,1) - 4);2*mat(i,2) + 2*mat(i,2)*(mat(i,1) - 2)^2 + 2];
    H=[12*(mat(i,1) - 2)^2 + 2*mat(i,2)^2, 4*(mat(i,1)-2)*mat(i,2);4*(mat(i,1)-2)*mat(i,2), 2*(mat(i,1) - 2)^2 + 2];
    p=H\J;
    mat(i+1,:) = mat(i,:)-p';
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end
end