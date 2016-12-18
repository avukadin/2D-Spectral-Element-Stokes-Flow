function [ polynomials] = nodal_coord(y1e,y2e,x1e,x2e,x_domain,y_domain,element,P,X,shift,C)

polynomials=zeros(size(x_domain,2),size(y_domain,2));

xt=@(zeta) .5*(x1e+x2e)+.5.*(x2e-x1e).*zeta;
Ldx1= (x_domain-xt(P(2))).*(x_domain-xt(P(3))).*(x_domain-xt(P(4))) ./( (xt(P(1))-xt(P(2))).*(xt(P(1))-xt(P(3))).*(xt(P(1))-xt(P(4))) );
Ldx2= (x_domain-xt(P(1))).*(x_domain-xt(P(3))).*(x_domain-xt(P(4))) ./( (xt(P(2))-xt(P(1))).*(xt(P(2))-xt(P(3))).*(xt(P(2))-xt(P(4))) );
Ldx3= (x_domain-xt(P(2))).*(x_domain-xt(P(1))).*(x_domain-xt(P(4))) ./( (xt(P(3))-xt(P(2))).*(xt(P(3))-xt(P(1))).*(xt(P(3))-xt(P(4))) );
Ldx4= (x_domain-xt(P(2))).*(x_domain-xt(P(3))).*(x_domain-xt(P(1))) ./( (xt(P(4))-xt(P(2))).*(xt(P(4))-xt(P(3))).*(xt(P(4))-xt(P(1))) );
Lx={Ldx1,Ldx2,Ldx3,Ldx4};

yt=@(eta) .5*(y1e+y2e)+.5.*(y2e-y1e).*eta;
Ldy1=(y_domain-yt(P(2))).*(y_domain-yt(P(3))).*(y_domain-yt(P(4))) ./( (yt(P(1))-yt(P(2))).*(yt(P(1))-yt(P(3))).*(yt(P(1))-yt(P(4))) );
Ldy2=(y_domain-yt(P(1))).*(y_domain-yt(P(3))).*(y_domain-yt(P(4))) ./( (yt(P(2))-yt(P(1))).*(yt(P(2))-yt(P(3))).*(yt(P(2))-yt(P(4))) );
Ldy3=(y_domain-yt(P(2))).*(y_domain-yt(P(1))).*(y_domain-yt(P(4))) ./( (yt(P(3))-yt(P(2))).*(yt(P(3))-yt(P(1))).*(yt(P(3))-yt(P(4))) );
Ldy4=(y_domain-yt(P(2))).*(y_domain-yt(P(3))).*(y_domain-yt(P(1))) ./( (yt(P(4))-yt(P(2))).*(yt(P(4))-yt(P(3))).*(yt(P(4))-yt(P(1))) );
Ly={Ldy1,Ldy2,Ldy3,Ldy4};

for q=1:size(x_domain,2)
    
    for r=1:size(y_domain,2)       
        node=1;
        
        for k=1:4
            
            for n=1:4   
                
                polynomials(r,q)=X(C(element,node)+shift)*Ly{n}(r)*Lx{k}(q)+polynomials(r,q);
                node=node+1;
                
            end
            
        end
        
    end

end

return;

end