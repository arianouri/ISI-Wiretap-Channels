function dCheck(MAT1,MAT2,intP1,intP2)
clc
    reshape(MAT1(intP1+1,intP2+1,1,:),1,2)
    reshape(MAT2(intP1+1,intP2+1,1,:),1,2)
end