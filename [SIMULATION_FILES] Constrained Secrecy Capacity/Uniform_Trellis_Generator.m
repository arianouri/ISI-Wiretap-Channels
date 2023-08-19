function [ jPR_Trellis ] = Uniform_Trellis_Generator(  Branch_Cardinal , Trellis_Extended_io)
jPR_Trellis=zeros(length(SpecialiZe(Trellis_Extended_io)),length(SpecialiZe(Trellis_Extended_io)));
 for iIndex=1:size(Trellis_Extended_io,1)
   jPR_Trellis(Trellis_Extended_io(iIndex,1),Trellis_Extended_io(iIndex,2))=1./length(SpecialiZe(Trellis_Extended_io(Trellis_Extended_io(:,1)==1,2)))./Branch_Cardinal;
 end
end

