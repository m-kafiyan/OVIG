function out = B(W0, v0, D)
mul=1;
for i=1:D
    mul = mul* gamma((v0+1-i)/2);
end
if  mul == Inf
    mul=1e100;
end
out = (det(W0)^(-v0/2)) * (((2^(v0*D/2))*(pi^(D*(D-1)/4))*mul)^(-1));
end