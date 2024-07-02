function autData = map_spectrum(lambda)
% MAP_SPECTRUM This function maps spectrum that is possibly complex to real
% spectrum. It can also be used to rearrange state variables.
% Given nonlinear programming is perfomred in real vector space. So
% the q here is arranged as p = [x1,real(x2),imag(x2),x3,...] provided x1
% is real and (x2,x3) is a complex conjugate pair, and x3 is real
% eigenvalues Lambda is sorted in descending order of real parts
% positive imaginary part is placed first in each complex pair

dim  = numel(lambda);
autData = struct();
autData.dim = dim;
% lambda = obj.System.spectrum.Lambda(modes);
realx  = find(abs(imag(lambda))<1e-6);
compx  = setdiff(1:dim, realx);
lamdReal = lambda(realx);
lamdComp = lambda(compx);
ncomp    = numel(compx);
assert(mod(ncomp,2)==0, 'Complex eigenvalues are not in pair');
lamdCompFirst  = lamdComp(1:2:end-1);
lamdCompSecond = lamdComp(2:2:end);
flag1 = abs(real(lamdCompFirst)-real(lamdCompSecond))<1e-6*abs(real(lamdCompFirst));
flag1 = all(flag1);
flag2 = abs(imag(lamdCompFirst)+imag(lamdCompSecond))<1e-6*abs(imag(lamdCompFirst));
flag2 = all(flag2);
flag3 = all(imag(lamdCompFirst)>0);
assert(flag1, 'Real parts of complex conjugate pair are not equal');
assert(flag2, 'Imaginary parts of complex conguate pair are not opposite');
assert(flag3, 'The imaginary part of the first entry in each pair is not positive');
lamd = zeros(dim,1);
lamd(realx) = real(lamdReal);
lamd(compx(1:2:end-1)) = real(lamdCompFirst);
lamd(compx(2:2:end))   = imag(lamdCompFirst);

lambda(compx(2:2:end)) = [];
idxReal = find(abs(imag(lambda))<1e-6);
idxComp = setdiff(1:numel(lambda), idxReal);

autData.realx   = realx;
autData.compx   = compx;
autData.lamd    = lamd;
autData.idxReal = idxReal;
autData.idxComp = idxComp;

end