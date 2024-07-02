%% Bernoulli beam with 1:3 internal resonance
% In this example, we consider a cantilever beam with a nonlinear support spring 
% at its free end. The linear part of the stiffness of the support spring is tuned 
% such that 1:3 internal resonance occurs between the first two bending modes.
% 
% 
% 
% We then extract the forced response curve for both periodic and quasi-periodic 
% response using SSM reduction. In particular, both two- and three-dimensional 
% invariant tori will be computed.
% Dynamical System Setup
% Numerical experiments show that a near 1:3 internal resonance occurs at $k_l=27$. 
% In the following computations, we set the number of beam elements to be 40. 
% The bifurcations observed here are persistent when the number of elements is 
% increased.

nElements = 40;
kLinear = 27;
kNonlinear = 60;
[M,C,K,fnl] = build_model(nElements,kLinear,kNonlinear);
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',10,'Nmax',10,'notation','multiindex','RayleighDamp',false)

%%
% Linear Modal analysis
%
[V,D,W] = DS.linear_spectral_analysis();

%% Add forcing
% Excitation of the form $\mathbf{f}=\omega_1^2\mathbf{M}\mathbf{\phi}_1\cos\Omega 
% t$ is applied such that only the first linear mode is activated if damping and 
% nonlinear internal forces are removed.
%

[vs,om_nat] = eigs(K,M,2,'smallestabs');
om_nat = sqrt(diag(om_nat));
f_0 = (om_nat(1))^2*M*vs(:,1);
epsilon = 0.002;
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas, epsilon);

%%
% SSM Computation
%
% *Choose master spectral subspace*
%

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')

resonant_modes = [1 2 3 4];
mFreq  = [1 3];
order  = 7; % Approximation order
outdof = [2*round(nElements/2)-1; 2*nElements-1];
n = length(M);

%% Continuation of equilibria
% SSM_isol2ep: continuation of equilibrium points
% Each equilibrium here corresponds to a periodic orbit of the full system

freqRange  = [15.30 15.95];
set(S.FRCOptions, 'SampStyle','cocoBD');                       % sampling style
set(S.FRCOptions, 'nCycle',5000, 'initialSolver', 'fsolve');   % solver for initial solution
set(S.contOptions, 'h_min',1e-3,'h_max',0.01,'NSV',10,'bi_direct',true,'PtMX',100);  % continuation setting
set(S.FRCOptions, 'coordinates', 'cartesian');                 % two representations

isolid = ['isol-',num2str(nElements),'-',num2str(order),'c'];
startep = tic;
FRC = S.SSM_isol2ep(isolid,resonant_modes,order,mFreq,'freq',freqRange,outdof);
timings.epFRC = toc(startep);

%% Continuation of Saddle Node bifurcations
% SSM_ep2SN: continuation of SN equilibrium points
% Continuation of saddle-node bifurcation of periodic orbits

set(S.contOptions, 'h_min',1e-3,'h_max',0.1);              % continuation setting
epsRange = [1e-4 1e-2];
bd    = coco_bd_read([isolid,'.ep']);
omega = coco_bd_col(bd,'om');
SNlab = coco_bd_labs(bd,'SN');
if ~isempty(SNlab)
    % find the lab with smallest omega
    SNidx = coco_bd_idxs(bd,'SN');
    omSN  = omega(SNidx);
    [~,id]= max(omSN);
    SNlab = SNlab(id);
    SNid = ['SN-',num2str(nElements),'-',num2str(order),'c'];
    S.SSM_ep2SN(SNid,isolid,SNlab,{freqRange,epsRange},outdof);
end

%% Continuation of Hopf Bifurcations
% SSM_ep2HB: continuation of HB equilibrium points
% Continuation of Hopf bifurcation of periodic orbits

HBlab = coco_bd_labs(bd,'HB');
% find the lab with smallest omega
HBidx = coco_bd_idxs(bd,'HB');
omHB  = omega(HBidx);
[~,idx] = sort(omHB);
HBlab1 = HBlab(idx(end-1));
HBid1 = ['HB1-',num2str(nElements),'-',num2str(order),'c'];
S.SSM_ep2HB(HBid1,isolid,HBlab1,{freqRange,epsRange},outdof);

%% Investigation of a specific Hopf Bifurcation
% Continuation of HB points with the other one as starting point
%

HBid2 = ['HB2-',num2str(nElements),'-',num2str(order),'c'];
S.SSM_ep2HB(HBid2,isolid,HBlab(idx(2)),{freqRange,epsRange},outdof);

%%
% SSM_HB2po: continuation of periodic orbits from HB point 

set(S.contOptions, 'h_max',0.3,'PtMX',100, 'bi_direct', false, 'NSV', 1,'NAdapt',5);                    % continuation setting
po1id = ['po1-',num2str(nElements),'-',num2str(order),'c'];
startpo = tic;
set(S.FRCOptions,'parSamps',15.75);
S.SSM_HB2po(po1id,isolid,HBlab1,'freq',freqRange,[outdof; outdof+n],'saveICs');
timings.po1FRC = toc(startpo);
%% 
% Continuation of periodic orbits from another HB point

HBlab2 = HBlab(idx(2));
po2id = ['po2-',num2str(nElements),'-',num2str(order),'c'];
set(S.contOptions, 'h_max',0.05, 'PtMX', 70, 'NSV', 2, 'bi_direct', false, 'NAdapt', 10);              % continuation setting
startpo = tic;
set(S.FRCOptions,'parSamps',[15.589, 15.5905]);
S.SSM_HB2po(po2id,isolid,HBlab2,'freq',freqRange,[outdof; outdof+n],'saveICs');
timings.po2FRC = toc(startpo);

%% Continuation of Torus Bifurcations
%  SSM_po2TR: continuation of TR bifurcation periodic orbits
% Continuation of quasi-periodic Hopf bifurcation of two-dimensional invairant 
% tori
%

bd    = coco_bd_read([po2id,'.po']);
TRlab = coco_bd_labs(bd,'TR');
assert(~isempty(TRlab), 'No TR periodic orbits are found');
set(S.contOptions, 'h_max',1, 'bi_direct', true, 'NAdapt', 0);              % continuation setting
TRid = ['TR-',num2str(nElements),'-',num2str(order),'c'];
S.SSM_po2TR(TRid,po2id,TRlab,{freqRange,epsRange},[outdof; outdof+n]);

%%
% SSM_TR2tor: continuation of tori from TR point
% Continuation of three-dimensional invariant tori
%


TRlab = 1; %
set(S.contOptions, 'h_max',100,'PtMX',50,'bi_direct',false,'NSV', 5);              % continuation setting
set(S.FRCOptions,"torNumSegs",15); %2*15+1=31 segments
torid = ['tor-',num2str(nElements),'-',num2str(order),'c'];
starttor = tic;
S.SSM_TR2tor(torid,TRid,TRlab,'freq',freqRange,[outdof; outdof+n],'saveICs');
timings.torFRC = toc(starttor);