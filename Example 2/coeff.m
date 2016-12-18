function [ coeff_vec ] = coeff( f,C,ne_y,ne_x,x1,y1,P,x_len,y_len,Ng)
%Returns the interpolant values of the function f.


coeff_vec=zeros(Ng,1);
element=1;

for k=1:ne_x
    
    x1e=x1+x_len*(k-1);
    x2e=x1+x_len*k;
    
    for p=1:ne_y
        
        y1e=y1+y_len*(p-1);
        y2e=y1+y_len*p;
        x=@(zeta) .5*(x1e+x2e)+.5*(x2e-x1e).*zeta;
        y=@(eta) .5*(y1e+y2e)+.5*(y2e-y1e).*eta;
        n=1;
        
        for i=1:4
            
            for j=1:4
                
                coeff_vec(C(element,n))= f(x(P(i)),y(P(j)));
                n=n+1;
                
            end
            
        end
        
        element=element+1;
        
    end
    
end

return;

end