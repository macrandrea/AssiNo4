maxit=10;
tol=1e-8;
mat=zeros(maxit+1,2);
mat(1,2)=1;
mat(1,1)=1;
j=1;
for i=1:maxit
    J=[4*(mat(i,1) - 2)^3 + mat(i,2)^2*(2*mat(i,1) - 4);2*mat(i,2) + 2*mat(i,2)*(mat(i,1) - 2)^2 + 2];
    H=[12*(mat(i,1) - 2)^2 + 2*mat(i,2)^2, 4*(mat(i,1)-2)*mat(i,2);4*(mat(i,1)-2)*mat(i,2), 2*(mat(i,1) - 2)^2 + 2];
    p=H\J;
    mat(i+1,:) = mat(i,:)-p';
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end
    j=j+1;
    %Ef(i)=(mat(i,1)-2)^4+(mat(i,1)-2)^2*mat(i,2)^2+(mat(i,2)+1)^2;
    hold on
    %plot([j-1 j], [mat(j-1,:) mat(j,:) ], 'r-');
    %fprintf('\n Valor= %8.3f  ',mat(i,:))
    %plot3(Ef(i),mat(i,1),Ef(i),mat(i,2), 'r-')
    %fplot (0)
    hold off
end