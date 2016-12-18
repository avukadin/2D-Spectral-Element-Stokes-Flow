clear variables;
%close all;
clc;
beep off;
figure('units','normalized','position',[.2 .2 .7 .7]);

%==========================
% Domain to be Discretised
%==========================
x1=0;
x2=2;
y1=0;
y2=2;
ne_x=10; % # of elements in x direction
ne_y=10; % # of elements in y direction
x_len=(x2-x1)/ne_x; % length of an element in x direction
y_len=(y2-y1)/ne_y; % length of an element in y direction
xres=.005; % x resolution
yres=xres*y_len/x_len; % y resolution
x_domain = x1:xres:x2;
y_domain = y1:yres:y2;
Ne=ne_x*ne_y; % Total Elements
J=x_len*y_len/4; %Jacobian

%=======================================
%Problem Setup
%=======================================
[C,boundary,Ng]=connectivity(ne_x,ne_y,Ne);% C: connectivity matrix || boundary: global boundary coordinates || Ng:# of global nodes
[P,Lu,L, L_dx, L_dy]=gll_lag(); %Gets GLL points and Lagrange basis polynomials and their derivatives
mu=1; %Dynamic Viscocity
rho=1; %Density
nu=mu/rho; %Kinematic Viscocity
step=1/10; %Time Step

%==================
%Initial Conditions
%==================
U_n=zeros(Ng*2,1);
u_init=@(x,y) 0;
v_init=@(x,y) 0;
u_coeff=coeff( u_init,C,ne_y,ne_x,x1,y1,P,x_len,y_len,Ng);
v_coeff=coeff( v_init,C,ne_y,ne_x,x1,y1,P,x_len,y_len,Ng);
for i=1:Ng
    U_n(i)=u_coeff(i);
    U_n(i+Ng)=v_coeff(i);
end

%=========================
%Gets all Global Matrices
%=========================
[gD_x,gD_y, gDu,gMu, gQ_x, gQ_y ,gR]=matrix_sk(L ,L_dx,L_dy,J,x_len,y_len,Ne,Ng,C);


frame=1;
for time=0:step:10
    b=zeros(Ng*2+Ne,1);
    G01=zeros(Ng,Ng);
    G02=zeros(Ne,Ne);
    G03=zeros(Ne,Ng);
    
    %=====================
    %Assmble Global System
    %=====================
    Global1=[mu.*((1/(nu*step)).*gMu+.5.*gDu)+.5*mu.*gR, G01, -gD_x'-gQ_x
        G01, mu.*((1/(nu*step)).*gMu+.5.*gDu)+.5*mu.*gR, -gD_y'-gQ_y
        -gD_x, -gD_y ,G02];
    Global2=[mu.*((1/(nu*step)).*gMu-.5.*gDu), G01
        G01,mu.*((1/(nu*step)).*gMu-.5.*gDu)
        G03, G03];
    b=Global2*U_n;
    
    %===================
    %Boundary Conditions
    %===================
    for j=1:Ng
        found=0;
        for k=1:size(boundary,2)
            
            if boundary(k)==j  
                found=1;
            end
        end
        
        if j<=3*ne_y+1
            if time/10<=1 %Boundary Velocity at Left Wall
                V=time/10;
                
            else
                V=1;
                
            end
            
        elseif j>=Ng-3*ne_y
            V=0;%Boundary Velocity at Left Wall
        
        else      
            V=0;
        end
        
        if found==1
            
            for i=1:2*Ng+Ne
                b(i)=b(i)-Global1(i,j)*0-Global1(i,Ng+j)*V;
                Global1(i,j)=0; Global1(i,Ng+j)=0;
                Global1(j,i)=0; Global1(Ng+j,i)=0;
            end    
            Global1(j,j)=1; Global1(j+Ng,j+Ng)=1;
            b(j)=0; b(Ng+j)=V;      
        end
    end
    
    
    %===========
    %Solve for X
    %===========
    X=Global1\(b);
    n=1;
    X(size(X,1))=1;
    pause(0.0001)
    cla;
    
    for i=1:ne_x
        %================================
        %Get x domain for current element
        %================================
        x1e=x1+(i-1)*x_len;
        x2e=x1+(i)*x_len;
        x_domain=x1e:xres:x2e;
        
        for j=1:ne_y
            %================================
            %Get y domain for current element
            %================================
            y1e=y1+(j-1)*y_len;
            y2e=y1+(j)*y_len;
            y_domain=y1e:yres:y2e;
            axis([0 2 0 2])
            
            velocity_u=nodal_coord(y1e,y2e,x1e,x2e,x_domain,y_domain,n,P,X,0,C);
            velocity_v=nodal_coord(y1e,y2e,x1e,x2e,x_domain,y_domain,n,P,X,Ng,C);
            
            hold on
            n=n+1;
            
            %==============
            %Plot Solution
            %==============
            scale=.2;
            [x_domain1,y_domain1] = meshgrid (x_domain,y_domain);
            str = sprintf('Example 3.1, Time = %f s' ,time);
            axis([0 2 0 2])
            title(str);
            xlabel('x');
            ylabel('y');
            caxis([0, 1]);
            %colorbar;
            colormap jet;
            %rotate3d on;
            %quiver(x_domain1,y_domain1,velocity_u.*scale,velocity_v.*scale,'Autoscale','off','color',[0 0 1]);
            speed=sqrt(velocity_v.^2+velocity_u.^2);
            surf(x_domain1,y_domain1,speed,'EdgeColor', 'none');
        end
    end
    
    time %Current Time Being Processed
    
    hold off
    
    
    mov(frame)=getframe(gcf,[0,0,1300,800]);%Creates Movie
    frame=frame+1;
    
    %=================================
    %Set Conditions for next iteration
    %=================================
    for i=1:2*Ng
        U_n(i)=X(i);
    end
    
end

%=================
%Export Movie File
%=================
mov(:,1)=[];
myVideo = VideoWriter('q33.avi');
myVideo.FrameRate = 10; 
myVideo.Quality = 100;    
open(myVideo);
writeVideo(myVideo, mov);
close(myVideo);