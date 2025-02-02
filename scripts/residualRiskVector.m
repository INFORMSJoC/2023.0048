function [ r ] = residualRiskVector( Omega, B, ki)


r = exp(-ki.*B);
r = repmat(r,size(Omega,1),1) .* Omega;
r = sum(r,2);


end

