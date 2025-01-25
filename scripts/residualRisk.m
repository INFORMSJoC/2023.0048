function [ r ] = residualRisk( mu, BS, ki)


r = exp(-repmat(ki,size(BS,1),1).*BS);
r = repmat(mu,size(BS,1),1) .* r;
r = sum(r,2);


end

