function dsdt = compsEqs(t,s,P,Bhat,Rhat,AHat,bQ,b)
% COMPSEQS This function computes the vector field of the compensated
% vector s.
%
% dsdt = compsEqs(t,s,P,Bhat,Rhat,LamHat,bQ,b)
%
% dsdt = (P*Bhat*R^(-1)*Bhat.')*s+bQ+b

Pt   = P(t); 
nq   = size(Bhat,1);
Pt   = reshape(Pt,[nq,nq]); % vector to matrix
dsdt = (Pt*Bhat*(Rhat\(Bhat.'))-AHat.')*s+bQ(t)+2*Pt*b(t);

end