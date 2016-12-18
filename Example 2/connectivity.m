function [ C,boundary,Ng] = connectivity( ne_x,ne_y,elements)

C=zeros(elements,16);%Connectivity Matrix
element=1;
start=1;
k=1;
boundary=zeros();

%=============================
%Construct Connectivity Matrix
%=============================

for n=1:ne_x
    
    for m=1:ne_y
        
        if m~=1
            
            start=C(element-1,4);
            
        elseif m==1 && n~=1
            
            start=C(element-ne_y,13);
            
        end
        
        for i=1:4
            
            C(element,i)= start+i-1;
            C(element,i+4)= start+4*(ne_y)-ne_y+1+i-1;
            C(element,i+8)= start+2*(4*(ne_y)-ne_y+1)+i-1;
            C(element,i+12)= start+3*(4*(ne_y)-ne_y+1)+i-1;
            
        end
        
        if m==ne_y
            
            boundary(size(boundary,2)+1)=C(element,4);
            boundary(size(boundary,2)+1)=C(element,8);
            boundary(size(boundary,2)+1)=C(element,12);
            boundary(size(boundary,2)+1)=C(element,16);
            k=k+4;
            
        end
        
        if m==1
            
            boundary(size(boundary,2)+1)=C(element,1);
            boundary(size(boundary,2)+1)=C(element,5);
            boundary(size(boundary,2)+1)=C(element,9);
            boundary(size(boundary,2)+1)=C(element,13);
            k=k+4;
            
        end
        
        if n==ne_x
            
            boundary(size(boundary,2)+1)=C(element,13);
            boundary(size(boundary,2)+1)=C(element,14);
            boundary(size(boundary,2)+1)=C(element,15);
            boundary(size(boundary,2)+1)=C(element,16);
            k=k+4;
            
        end
        
        if n==1
            
            boundary(size(boundary,2)+1)=C(element,1);
            boundary(size(boundary,2)+1)=C(element,2);
            boundary(size(boundary,2)+1)=C(element,3);
            boundary(size(boundary,2)+1)=C(element,4);
            k=k+4;
            
        end
        
        element=element+1;
        
    end
    
end

boundary(1)=boundary(2);
boundary=unique(sort(boundary));
Ng=C(size(C,1),size(C,2));% number of Global elements

return;

end