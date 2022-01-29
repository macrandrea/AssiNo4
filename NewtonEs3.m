clc,clear
maxit=50;
tol=1e-8;
mat=zeros(maxit+1,2);
alpha=ones(maxit,1);
x0=[8,0.8];
mat(1,:)=x0;
delta=0.25;
epsilon=0.75;
for i=1:maxit
    H=[2*(1 - mat(i,2)^2)^2 + 2*(1 - mat(i,2)^3)^2, 2* (1 - mat(i,2)) - 4 *mat(i,1)* mat(i,2)* (1 - mat(i,2)^2) - 6* mat(i,1)* mat(i,2)^2* (1 - mat(i,2)^3) + 4 *mat(i,2) *(2.25 - mat(i,1)* (1 - mat(i,2)^2)) + 6 *mat(i,2)^2* (2.625 - mat(i,1) *(1 - mat(i,2)^3)); 2 *(1 - mat(i,2)) - 4* mat(i,1)* mat(i,2) *(1 - mat(i,2)^2) - 6 *mat(i,1) *mat(i,2)^2 *(1 - mat(i,2)^3) + 4 *mat(i,2) *(2.25 - mat(i,1)* (1 - mat(i,2)^2)) + 6 *mat(i,2)^2 *(2.625 - mat(i,1) *(1 - mat(i,2)^3)), -2* mat(i,1) + 8* mat(i,1)^2 *mat(i,2)^2 + 18 *mat(i,1)^2 *mat(i,2)^4 + 4* mat(i,1) *(2.25 - mat(i,1) *(1 - mat(i,2)^2)) + 12 *mat(i,1) *mat(i,2)* (2.625 - mat(i,1) *(1 - mat(i,2)^3))];
    J=[2 *(mat(i,2)^3 - 1)* (2.625 - mat(i,1)* (1 - mat(i,2)^3)) + 2* (mat(i,2)^2 - 1)* (2.25 - mat(i,1)* (1 - mat(i,2)^2)) - (1 - mat(i,2))^2 ; 4* mat(i,1)* mat(i,2)* (2.25 - mat(i,1)* (1 - mat(i,2)^2)) + 6* mat(i,1)* mat(i,2)^2* (2.625 - mat(i,1)* (1 - mat(i,2)^3)) + 2* mat(i,1)* (1 - mat(i,2))];
    p=H\J;
    mat(i+1,:) = mat(i,:)-alpha(i)*p';
    EfPO=(1.5-mat(i+1,1)*(1-mat(i+1,2))^2)+(2.25-mat(i+1,1)*(1-mat(i+1,2)^2))^2+(2.625-mat(i+1,1)*(1-mat(i+1,2)^3))^2;
    Ef  =(1.5-mat(i,1)*(1-mat(i,2))^2)+(2.25-mat(i,1)*(1-mat(i,2)^2))^2+(2.625-mat(i,1)*(1-mat(i,2)^3))^2;
    if (EfPO > Ef+ delta*alpha*sum(J'*p))
          alpha(i+1)=epsilon*alpha(i); 
    else
          alpha(i+1)=alpha(i);
        %if (alpha(i)<0.19)
        %   
        %end
    end
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end   
end

% f= @(x,y) (1.5-x*(1-y)^2)+(2.25-x*(1-y^2))^2+(2.625-x*(1-y^3))^2;
% hes=hessian(f);
