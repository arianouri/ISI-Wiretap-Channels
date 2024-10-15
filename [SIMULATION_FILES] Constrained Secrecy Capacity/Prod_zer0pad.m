function [Ret1] = Prod_zer0pad(A,B)
    auxA=[A zeros(1,length(B)-length(A))];
    auxB=[B zeros(1,length(A)-length(B))];
    Ret1=auxA.*auxB;
end