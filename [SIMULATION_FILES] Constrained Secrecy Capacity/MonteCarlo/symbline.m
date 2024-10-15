function SymbLine=symbline(Length,Symbol)
SymbLine='';
for Index=1:Length
    SymbLine=[SymbLine Symbol];
end
disp(SymbLine)