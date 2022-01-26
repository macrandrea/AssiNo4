clc,clear
maxit=10;
tol=1e-8;
x0=[-1,3,3,0]';
b=[5.04, -59.4, 146.4, -96.6]';
A=[ 0.16  , -1.2 , 2.4 , -1.4;
    -1.2 , 12.0 , -27 , 16.8;
     2.4 , -27.0, 64.8, -42.0;
     -1.4, 16.8 ,-42.0, 28.0];
mat=zeros(maxit,4);
mat(1,:)=x0;
r=zeros(maxit,4);
p=r;
alpha=zeros(maxit+1);
beta=alpha;
r(1,:)=b-A*x0;
p(1,:)=r(1,:);
for i=1:(maxit+1)
    alpha(i)=(norm(r(i,:))^2)/(p(i,:)*A*p(i,:)');
    mat(i+1,:)=mat(i,:)+alpha(i)*p(i,:);
    r(i+1,:)=b-A*mat(i+1,:)';
    beta(i)=norm(r(i+1,:))^2/norm(r(i,:))^2;
    p(i+1,:)=r(i+1,:)+beta(i)*p(i,:);
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))))
        break;
    end
end