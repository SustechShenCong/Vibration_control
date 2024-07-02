function  dPdt = RiccatiEqs(t,P,A,B,Q,R)
% RICCATIEQS: This function returns the vector field for Riccati
% differential equations. Here X is a column vector obtained from the
% vectorization of matrix
%
% dPdt = RiccatiEqs(t,P,A,B,Q,R)
%
% dPdt = -P*A-A.'*P-Q+P*B*R^(-1)*B.'*P

P    = reshape(P,size(A)); % vector to matrix
dPdt = -P*A-A.'*P-Q+P*B*(R\(B.'*P));
dPdt = dPdt(:);            % matrix to vector      

end