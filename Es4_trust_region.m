maxit=20;
tol=1e-8;
eta=0.01;
mat=zeros(maxit+1,2);
alpha=ones(maxit,1);
f=zeros(maxit+1,1);
mat(1,:)=[0,0];
x=mat(1,1);
y=mat(1,2);
f(1)=x^4+x*y+(1+y)^2;
flag=0;
for i=1:maxit
    H=[12*x^2, 1;
       1 , 2];
    J=[4*x^3+y; x+2*y+2];
    [Q,L]=eig(H);
    mineig=min(diag(L));
    u=0;
    if mineig<=0
        u=1.01*abs(mineig)+0.01;
    end
    q=Q'*J;
    rho=0;
    while rho<eta
        squarenormp=abs(q')*abs((diag(L)+u).^-1);
        if flag==0
            flag=1;
            delta=squarenormp;
        end
        while squarenormp>delta^2
            if u==0
                u=0.01;
            else
                u=2*u;
            end
            squarenormp=abs(q')*abs((diag(L)+u).^-1);
        end
        p=Q*diag((diag(L)+u).^-1)*Q'*J;
        mat(i+1,:) = mat(i,:)-p';
        x=mat(i+1,1);
        y=mat(i+1,2);
        f(i+1)=x^4+x*y+(1+y)^2;
        rho=(f(i)-f(i+1))/(3/2*p'*J);
        if rho<0.25
            delta=delta/4;
        elseif rho>0.75
            delta=delta*2;
        end
    end
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end
end
