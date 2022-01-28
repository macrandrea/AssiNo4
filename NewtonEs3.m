clc,clear
maxit=100;
tol=1e-8;
mat=zeros(maxit+1,2);
alpha=ones(maxit,1);
x0=[8,0.2];
mat(1,:)=x0;
delta=0.5;
epsilon=0.5;
for i=1:maxit
    J=[4*(mat(i,1) - 2)^3 + mat(i,2)^2*(2*mat(i,1) - 4);2*mat(i,2) + 2*mat(i,2)*(mat(i,1) - 2)^2 + 2];
    H=[12*(mat(i,1) - 2)^2 + 2*mat(i,2)^2, 4*(mat(i,1)-2)*mat(i,2);4*(mat(i,1)-2)*mat(i,2), 2*(mat(i,1) - 2)^2 + 2];
    p=H\J;
    mat(i+1,:) = mat(i,:)-alpha(i)*p';
    EfPO=(1.5-mat(i+1,1)*(1-mat(i+1,2))^2)+(2.25-mat(i+1,1)*(1-mat(i+1,2)^2))^2+(2.625-mat(i+1,1)*(1-mat(i+1,2)^3))^2;
    Ef  =(1.5-mat(i,1)*(1-mat(i,2))^2)+(2.25-mat(i,1)*(1-mat(i,2)^2))^2+(2.625-mat(i,1)*(1-mat(i,2)^3))^2;
    if(EfPO< Ef+delta*alpha(i)*sum(J'*p))
        alpha(i+1)=alpha(i)*epsilon;
        if (alpha(i)<0.1)
            alpha(i)=alpha(i-1);
        end
    end
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end
end