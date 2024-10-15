function [ SpeliZ ] = SpecialiZe( inpMatrix )
    SpeliZ=zeros(size(inpMatrix));
    sElements=1;
    sIndex=size(inpMatrix);
    for sIndexSize=1:size(sIndex,2)
        sElements=sElements*sIndex(sIndexSize);
    end
    SpIndex=0;
    for sIndex=1:sElements
    if sIndex==1
        SpIndex=SpIndex+1;
        SpeliZ(SpIndex)=inpMatrix(sIndex);
        sIndex=sIndex+1;
    end
    if ~isscalar(SpeliZ(SpeliZ==inpMatrix(sIndex)))
        SpIndex=SpIndex+1;
        SpeliZ(SpIndex)=inpMatrix(sIndex);
    end
    end
end

