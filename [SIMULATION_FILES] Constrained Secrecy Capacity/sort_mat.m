function Atld = sort_mat(A,colindex,aprc)
% aprc=1 -> Increasing / aprc=0 -> Decreasing
% Select r/c argument to 2 for column arrangement
    rcArg=1;
    Atld=A;
    for i=2:size(Atld,rcArg)
       X=Atld(i,:);
       x=Atld(i,colindex);
       j=i-1;
       if aprc==0
           while((j>=1)&&(abs(Atld(j,colindex))<abs(x)))
            Atld(j+1,:)=Atld(j,:);
            j=j-1;
           end
       elseif aprc==1
           while((j>=1)&&(Atld(j,colindex)>x))
            Atld(j+1,:)=Atld(j,:);
            j=j-1;
           end
       end
       Atld(j+1,:)=X;
    end
end

