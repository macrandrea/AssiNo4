%tic
%clear
%clc
%x=1;
%y=1;
%J=[4*(x - 2)^3 + y^2*(2*x - 4);2*y + 2*y*(x - 2)^2 + 2];
%H=[12*(x - 2)^2 + 2*y^2, 4*(x-2)*y;4*(x-2)*y, 2*(x - 2)^2 + 2];
%n=10;
%mat=zeros(n,2);
%mat(1,2)=y;
%mat(1,1)=x;
%j=0;
%while(norm(J)>0.01)
% for i=2:n
%        J=[4*(mat(i,1) - 2)^3 + mat(i,2)^2*(2*mat(i,1) - 4);2*mat(i,2) + 2*mat(i,2)*(mat(i,1) - 2)^2 + 2];
%        H=[12*(mat(i,1) - 2)^2 + 2*mat(i,2)^2, 4*(mat(i,1)-2)*mat(i,2);4*(mat(i,1)-2)*mat(i,2), 2*(mat(i,1) - 2)^2 + 2];
%        p=H\J;
%        mat(i,:) = mat(i,:) - p';
%        if (norm(mat(i,:)-mat(i-1,:))<0.01*(1+norm(mat(i,:))))
%            break;
%        end
% end
% j=j+1;
%end
%toc
x=[2,-1];
delta=zeros(10,1);
G=delta;
H=delta;
Ef=delta;
for i=1:10
delta(i)=norm(mat(i,:)-x);
J=[4*(mat(i,1) - 2)^3 + mat(i,2)^2*(2*mat(i,1) - 4);2*mat(i,2) + 2*mat(i,2)*(mat(i,1) - 2)^2 + 2];
H=[12*(mat(i,1) - 2)^2 + 2*mat(i,2)^2, 4*(mat(i,1)-2)*mat(i,2);4*(mat(i,1)-2)*mat(i,2), 2*(mat(i,1) - 2)^2 + 2];
G(i)=-J'*H^(-1)*J;
Ef(i)=(mat(i,1)-2)^4+(mat(i,1)-2)^2*mat(i,2)^2+(mat(i,2)+1)^2;
end
f = @(x,y) (x-2)^4+(x-2)^2*y^2+(y+1)^2;
fsurf(f,[-4 4 -4 4],"MeshDensity",100)
%fcontour(f,[-4 4 -4 4],"LevelStep",i)
 hold on
 plot([2,2],[-1,-1],'o')
 plot(mat(1:9,1),mat(1:9,2),'h-')
 %plot(mat,'red')
 %plot(Ef,'red')
 hold off