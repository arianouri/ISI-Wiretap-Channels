% for xIndex=1:8
%     for yIndex=2*xIndex-1:2*xIndex
%         
%         [xIndex mod(yIndex,8)]
%     end
% end

function [ jPR_Trellis ] = jSourceChannel( PR_Trellis , TRANS )
jPR_Trellis=zeros(size(TRANS));
for rIndex=1:size(TRANS,1)
    for cIndex=2*rIndex-1:2*rIndex
        I=mod(rIndex,size(PR_Trellis,1));
            I(I==0)=size(PR_Trellis,1);
        J=mod(cIndex,size(PR_Trellis,1));
            J(J==0)=size(PR_Trellis,2);
            ccIndex=mod(cIndex,size(TRANS,1));
            ccIndex(ccIndex==0)=size(TRANS,1);
        jPR_Trellis(rIndex,ccIndex)=PR_Trellis(I,J);
    end
end
end

