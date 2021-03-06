clc,clear
maxit=50;
tol=1e-8;
mat=zeros(maxit+1,2);
alpha=ones(maxit,1);
x0=[8,0.2];
mat(1,:)=x0;
delta=0.25;
epsilon=0.75;
eta=0.02;
mu=0.5;
EfZero=(1.5-0*(1-0)^2)+(2.25-0*(1-0^2))^2+(2.625-0*(1-0^3))^2;
HZero=[2*(1 - 0^2)^2 + 2*(1 - 0^3)^2, 2* (1 - 0) - 4 *0* 0* (1 - 0^2) - 6* 0* 0^2* (1 - 0^3) + 4 *0 *(2.25 - 0* (1 - 0^2)) + 6 *0^2* (2.625 - 0 *(1 - 0^3)); 2 *(1 - 0) - 4* 0* 0 *(1 - 0^2) - 6 *0 *0^2 *(1 - 0^3) + 4 *0 *(2.25 - 0* (1 - 0^2)) + 6 *0^2 *(2.625 - 0 *(1 - 0^3)), -2* 0 + 8* 0^2 *0^2 + 18 *0^2 *0^4 + 4* 0 *(2.25 - 0 *(1 - 0^2)) + 12 *0 *0* (2.625 - 0 *(1 - 0^3))];
JZero=[2 *(0^3 - 1)* (2.625 - 0* (1 - 0^3)) + 2* (0^2 - 1)* (2.25 - 0* (1 - 0^2)) - (1 - 0)^2 ; 4* 0* 0* (2.25 - 0* (1 - 0^2)) + 6* 0* 0^2* (2.625 - 0* (1 - 0^3)) + 2* 0* (1 - 0)];
em=zeros(maxit,1);
for i=1:maxit
    H=[2*(1 - mat(i,2)^2)^2 + 2*(1 - mat(i,2)^3)^2, 2* (1 - mat(i,2)) - 4 *mat(i,1)* mat(i,2)* (1 - mat(i,2)^2) - 6* mat(i,1)* mat(i,2)^2* (1 - mat(i,2)^3) + 4 *mat(i,2) *(2.25 - mat(i,1)* (1 - mat(i,2)^2)) + 6 *mat(i,2)^2* (2.625 - mat(i,1) *(1 - mat(i,2)^3)); 2 *(1 - mat(i,2)) - 4* mat(i,1)* mat(i,2) *(1 - mat(i,2)^2) - 6 *mat(i,1) *mat(i,2)^2 *(1 - mat(i,2)^3) + 4 *mat(i,2) *(2.25 - mat(i,1)* (1 - mat(i,2)^2)) + 6 *mat(i,2)^2 *(2.625 - mat(i,1) *(1 - mat(i,2)^3)), -2* mat(i,1) + 8* mat(i,1)^2 *mat(i,2)^2 + 18 *mat(i,1)^2 *mat(i,2)^4 + 4* mat(i,1) *(2.25 - mat(i,1) *(1 - mat(i,2)^2)) + 12 *mat(i,1) *mat(i,2)* (2.625 - mat(i,1) *(1 - mat(i,2)^3))];
    J=[2 *(mat(i,2)^3 - 1)* (2.625 - mat(i,1)* (1 - mat(i,2)^3)) + 2* (mat(i,2)^2 - 1)* (2.25 - mat(i,1)* (1 - mat(i,2)^2)) - (1 - mat(i,2))^2 ; 4* mat(i,1)* mat(i,2)* (2.25 - mat(i,1)* (1 - mat(i,2)^2)) + 6* mat(i,1)* mat(i,2)^2* (2.625 - mat(i,1)* (1 - mat(i,2)^3)) + 2* mat(i,1)* (1 - mat(i,2))];
    EfPO=(1.5-mat(i+1,1)*(1-mat(i+1,2))^2)+(2.25-mat(i+1,1)*(1-mat(i+1,2)^2))^2+(2.625-mat(i+1,1)*(1-mat(i+1,2)^3))^2;
    Ef  =(1.5-mat(i,1)*(1-mat(i,2))^2)+(2.25-mat(i,1)*(1-mat(i,2)^2))^2+(2.625-mat(i,1)*(1-mat(i,2)^3))^2;
    %p=H\J;
    
    p=-inv(H+mu*eye(2))*J;deltaTR=norm(p);
    emZero=EfZero+p'*JZero+0.5*p'*HZero*p;
    em(i)=Ef+p'*J+0.5*p'*H*p;
    rho=(Ef-EfPO)/(emZero-em(i));
    if ((0<rho)&&(rho<0.25))
        deltaTR=0.25*deltaTR;
    elseif (rho>0.75)
        deltaTR=2*deltaTR;
    elseif ((0.25<rho)&&(rho<0.75))
        
    end
    if (rho>eta)
    mat(i+1,:) = mat(i,:)+alpha(i)*p';
    else
        mat(i+1,:) = mat(i,:);
    end
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
