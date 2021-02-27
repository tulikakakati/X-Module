function score=ansim(g1,g2)

n=size(g1,2);
tempn=n;


score=0;




lmean=(g1(1)+g1(2)+g2(1)+g2(2))/4;
%disp('first column')
%disp(lmean);
 temp=2*max(abs(lmean-g1(1)),abs(lmean-g2(1)));
        
        if(temp~=0)
            score=score+(abs(g1(1)-g2(1))/temp);
        else
            tempn=tempn-1;
        end





for i=2:n-1
    
        lmean=(g1(i-1)+g1(i)+g1(i+1)+g2(i-1)+g2(i)+g2(i+1))/6;
    %disp(lmean);
        temp=2*max(abs(lmean-g1(i)),abs(lmean-g2(i)));
        
        if(temp~=0)
            score=score+(abs(g1(i)-g2(i))/temp);
        else
            tempn=tempn-1;
        end
       
    
end



lmean=(g1(n-1)+g1(n)+g2(n-1)+g2(n))/4;
%display('last column');
%disp(lmean);
 temp=2*max(abs(lmean-g1(n)),abs(lmean-g2(n)));
        
        if(temp~=0)
            score=score+(abs(g1(n)-g2(n))/temp);
        else
            tempn=tempn-1;
        end







score=1-(score/tempn);


end