maxit=1000;
tol=1e-8;
mat=zeros(maxit+1,4);
mat(1,:)=[-1,3,3,0]';
b=[5.04, -59.4, 146.4, -96.6]';
H=[ 0.16  , -1.2 , 2.4 , -1.4;
    -1.2 , 12.0 , -27 , 16.8;
     2.4 , -27.0, 64.8, -42.0;
     -1.4, 16.8 ,-42.0, 28.0];
for i=1:maxit
    J=b+H*mat(i,:)';
    p=H\J;
    mat(i+1,:) = mat(i,:)-p';
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end
end
