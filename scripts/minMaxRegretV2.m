

function [ MMR, B , c2, RR_star,binding,regretFunc] = minMaxRegretV2( MUs, BT, K)
% tic
Omega = generateScenarios(MUs);
numberOfDiseases = size(MUs,1);

Psi = size(Omega,1); % number of scenarios

RR_star = zeros(Psi,1);

c2 = true;
% tic
for i = 1:Psi
%     Omega(i,:)
    B_star = BSi_star_expo(Omega(i,:)',BT, K);
    if (sum(B_star<1e-7)>0)
%         B_star
        c2 = false;
    end
    RR_star(i) = residualRisk(Omega(i,:),B_star,K);
end
% toc
% RR_star = RR_star([1 2 4 5 6 7 8])
% Omega = Omega([1 2 4 5 6 7 8],:)
% RR_star
op = optimoptions('fminimax','Display','off','TolX',1e-11,'TolFun',1e-11,'TolCon',1e-11);
regr = @(x) Regret(x',RR_star,Omega,K);
x0 = ones(numberOfDiseases,1);
% x0 = zeros(numberOfDiseases,1);
% x0(1) = BT;
% if c2 == false
%     disp('Condition C2 not satisfied!');
% else
%     disp('Condition C2 IS satisfied!');
% end
% try
[B, MMR,~,~,~,lam] = fminimax(regr,x0,[],[],ones(1,numberOfDiseases),BT,zeros(numberOfDiseases,1),[],[],op);
% MMR = max(MMR);
B = B';
% lam
% lam.lower
% lam.upper
% catch
%     MUs
%     RR_star
%     Omega
%     K
%     regr
% end
binding = Omega(abs(MMR- max(MMR))<1e-7,:);
% toc
Regretmax = max(MMR);
regretFunc = regr;
end