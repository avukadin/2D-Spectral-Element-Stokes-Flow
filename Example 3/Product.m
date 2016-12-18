function [ prod ] = Product( f1,f2 )

prod=0;
w=[1/6,5/6,5/6,1/6];
P=[-1 ,-sqrt(1/5) , sqrt(1/5) , 1];

for i=1:4
    
    for j=1:4
        
        prod=prod+w(i)*w(j)*f1(P(i),P(j))*f2(P(i),P(j));
        
    end
    
end

end