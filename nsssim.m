function score=nsssim(g1,g2)

    n=size(g1,2);

    p1=zeros(1,n-1);
    p2=zeros(1,n-1);
    
    for i=2:n
        
        p1(1,i-1)=(g1(i)-g1(i-1))/(g1(2)-g1(1));
        p2(1,i-1)=(g2(i)-g2(i-1))/(g2(2)-g2(1));
        
    end
    
    score=ansim(p1(1,:),p2(1,:));
%disp(p1);
%disp(p2);
end