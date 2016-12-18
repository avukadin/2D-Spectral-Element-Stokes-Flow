function [ gD_x,gD_y, gDu,gMu, gQ_x, gQ_y ,gR] = matrix_sk( L ,L_dx,L_dy,J,x_len,y_len,Ne,Ng,C)
%Constructs the mass, diffusion and advection matrix
w=[1/6,5/6,5/6,1/6];
P=[-1 ,-sqrt(1/5) , sqrt(1/5) , 1];

D_x=zeros(16,1);
D_y=zeros(16,1);
Q_x=zeros(16,1);
Q_y=zeros(16,1);
gQ_x=zeros(Ng,Ne);
gQ_y=zeros(Ng,Ne);
R=zeros(16,16);
gD_x=zeros(Ne,Ng);
gD_y=zeros(Ne,Ng);
temp=@(x,y) 1;
gDu=zeros(Ng,Ng);
Du=zeros(16,16);
Mu=zeros(16,16);
gMu=zeros(Ng,Ng);
gR=zeros(Ng,Ng);
norm=1;
for i=1:16
    
    for j=1:16
        
        R(i,j)=norm.*(y_len/2).*(w(1)*w(1).*L{i}(P(1),P(1)).*(L_dx{j}(P(1),P(1)))+w(1)*w(2).*L{i}(P(1),P(2)).*(L_dx{j}(P(1),P(2)))+w(1)*w(3).*L{i}(P(1),P(3)).*(L_dx{j}(P(1),P(3)))+w(1)*w(4).*L{i}(P(1),P(4)).*(L_dx{j}(P(1),P(4)))...
            -norm.*(w(4)*w(1).*L{i}(P(4),P(1)).*(L_dx{j}(P(4),P(1)))+w(4)*w(2).*L{i}(P(4),P(2)).*(L_dx{j}(P(4),P(2)))+w(4)*w(3).*L{i}(P(4),P(3)).*(L_dx{j}(P(4),P(3)))+w(4)*w(4).*L{i}(P(4),P(4)).*(L_dx{j}(P(4),P(4)))))...
            +norm.*(x_len/2).*((w(1)*w(1).*L{i}(P(1),P(1)).*(L_dy{j}(P(1),P(1)))+w(1)*w(2).*L{i}(P(2),P(1)).*(L_dy{j}(P(2),P(1)))+w(1)*w(3).*L{i}(P(3),P(1)).*(L_dy{j}(P(3),P(1)))+w(1)*w(4).*L{i}(P(4),P(1)).*(L_dy{j}(P(4),P(1))))...
            -norm.*(w(4)*w(1).*L{i}(P(1),P(4)).*(L_dy{j}(P(1),P(4)))+w(4)*w(2).*L{i}(P(2),P(4)).*(L_dy{j}(P(2),P(4)))+w(4)*w(3).*L{i}(P(3),P(4)).*(L_dy{j}(P(3),P(4)))+w(4)*w(4).*L{i}(P(4),P(4)).*(L_dy{j}(P(4),P(4)))));  
    end  
end

for i=1:16
    D_x(i)=.5*y_len*Product(temp,L_dx{i});
    D_y(i)=.5*x_len*Product(temp,L_dy{i});   
    
    if i<=4
        Q_x(i)=norm.*(J*(w(1)*w(1)*L{i}(P(1),P(1)))+J*(w(1)*w(2)*L{i}(P(1),P(2)))+J*(w(1)*w(3)*L{i}(P(1),P(3)))+J*(w(1)*w(4)*L{i}(P(1),P(4))));  
    end
    
    if i>=13
        Q_x(i)=-norm.*(J*(w(4)*w(1)*L{i}(P(4),P(1)))+J*(w(4)*w(2)*L{i}(P(4),P(2)))+J*(w(4)*w(3)*L{i}(P(4),P(3)))+J*(w(4)*w(4)*L{i}(P(4),P(4))));
    end
    
    if i==4 || i==8 || i==12 || i==16
        Q_y(i)=-norm.*(J*(w(4)*w(1)*L{i}(P(1),P(4)))+J*(w(4)*w(2)*L{i}(P(2),P(4)))+J*(w(4)*w(3)*L{i}(P(3),P(4)))+J*(w(4)*w(4)*L{i}(P(4),P(4))));
    end
    
    if i==1 || i==5 || i==9 || i==13
        Q_y(i)=norm.*(J*(w(1)*w(1)*L{i}(P(1),P(1)))+J*(w(2)*w(1)*L{i}(P(2),P(1)))+J*(w(3)*w(1)*L{i}(P(3),P(1)))+J*(w(4)*w(1)*L{i}(P(4),P(1))));
    end
    
    for j=1:16
        Du(i,j)= (y_len/x_len)*Product(L_dx{i},L_dx{j})+ (x_len/y_len)*Product(L_dy{i},L_dy{j});
        Mu(i,j)=J*Product(L{i},L{j});  
    end
end

for i=1:Ne
    for j=1:16    
        gD_x(i,C(i,j))=D_x(j);
        gD_y(i,C(i,j))=D_y(j);
        gQ_x(C(i,j),i)=Q_x(j);%+gQ_x(C(i,j),i);
        gQ_y(C(i,j),i)=Q_y(j);%+gQ_y(C(i,j),i);  
        for k=1:16
            gDu(C(i,j),C(i,k))= Du(j,k) + gDu(C(i,j),C(i,k));
            gMu(C(i,j),C(i,k))= Mu(j,k)+gMu(C(i,j),C(i,k));
            gR(C(i,j),C(i,k))=R(j,k)+gR(C(i,j),C(i,k));
            
        end
        
    end
    
end

return;

end