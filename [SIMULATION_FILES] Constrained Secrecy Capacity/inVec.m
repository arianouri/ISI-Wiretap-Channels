function [ SaroTah ] = inVec( Vector )
SaroTah = zeros(size(Vector));
for i =1: size(Vector,2)
    SaroTah(i) = Vector(size(Vector,2)+1-i);
end;clear i
end

