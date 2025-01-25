function [ Omega] = generateScenarios( MUs )

n = size(MUs,1);

Omega = zeros(2^n,n);

for i = 1:n
    for j = 1:2^n
        index = rem(floor((j-1)/(2^n/2^i)),2)+1;
        Omega(j,i) = MUs(i,index);
    end
end


end
