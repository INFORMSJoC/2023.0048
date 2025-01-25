function [ R ] = Regret( B, R_star, Omega, K )

R = residualRiskVector( Omega, B, K)-R_star;

end

