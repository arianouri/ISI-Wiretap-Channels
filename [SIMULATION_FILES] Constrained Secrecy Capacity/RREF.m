function [X,REF,RREF]=RREF(Arg)
    REF=ones(size(Arg,1),size(Arg,2)+1)./0;
    ind=0;
    Atild=[Arg, [1:size(Arg,1)]'];
%     BROKER=1;
%     while((size(Atild,1)~=0)&&(BROKER))
    while size(Atild,1)~=0
      ind=ind+1;
      if((sum(Atild(1:end-1))>-10^(-8)&&(sum(Atild(1:end-1))<10^(-8))&&(size(Atild,1)==1)))
        REF(ind,:)=Atild(1,:);
        break;
      end
      [ArgCanbeCube,Atild]=pivot_fw_step(Atild);
      if size(Atild)==size(ArgCanbeCube)
%         BROKER=0;
        REF(ind:end,:)=ArgCanbeCube;
        break;
      end
      REF(ind,:)=ArgCanbeCube(1,:);
    end
    RREF = pivot_bw(REF);
%     RREF = sort_mat(RREF,size(RREF,2),1);
    X=RREF(:,size(RREF,2)-1);
end