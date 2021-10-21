function [seq]=equidistributed_generator(N,S_Cardinal)
%     N: length of the sequence
%     S_Cardinal: dimention;
    prime_vec=nthprime(1:S_Cardinal);
    seq=zeros(N,S_Cardinal);
    for inDex=1:N
        seq(inDex,:)=inDex.*sqrt(prime_vec);
        seq=mod(seq,1);
    end
end