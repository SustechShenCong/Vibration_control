function traj = LQR_closed_loop_BT(DS,z0,proj,tspan,autData,Wmap,masterModes,cont)
% LQR_closed_loops This routune presents LQR optimal control via SSM-based
% model reduction. Here we add closed-loop control strategy.

%% setup of system and control
Q    = cont.Q;
Rhat = cont.Rhat;
Mhat = cont.Mhat;
ep   = DS.epsilon;
n    = size(DS.M,1);

t0   = tspan(1); t1 = tspan(end);
if isempty(DS.E)
    Fext = @(t) zeros(2*n,1);
else
    Fext = @(t) [DS.E(t); zeros(n,1)];
end
%% setup initial conditions
p0   = get_initial_p0(DS,masterModes,z0,proj,autData,Wmap);
X0   = (z0-reduced_to_full(p0,Wmap,[],0))/ep;
%% Balanced Truncation
timeBT = tic;
Binv   = DS.B\eye(DS.N);
Atil   = Binv*DS.A;
Bext = [DS.D; zeros(size(DS.D))];
Btil   = Binv*Bext;
% Ctil   = eye(size(Atil,1));
% LX = lyapchol(Atil,Btil);
% LY = lyapchol(Atil',Ctil');
% [U, SIGMA, Z] = svd(LX*LY');
load("U.mat","U")
load("Z.mat","Z")
load("SIGMA.mat","SIGMA")
load("LX.mat","LX")
load("LY.mat","LY")
% disp("[Σ] = ")
% disp(SIGMA)
% tol_percentage = 0.6;
tol_percentage = 0.4;
k = 1;
% N = size(DS.A,1);
sigma_sum = SIGMA(1,1);
Sumation_SIGMA = sum(diag(SIGMA));
% 按相对误差的百分比，超过 tol_percentage % 认为误差会较大，在此处截断
while sigma_sum/Sumation_SIGMA < (1-tol_percentage)
    k = k + 1;
    sigma_sum = sigma_sum + SIGMA(k,k);
end

truncated_dim = k-1; % 截断维度，应由误差限决定
disp("truncated dimmension = ")
disp(truncated_dim)
fprintf('Time for Balanced truncation is %d\n',toc(timeBT));

sigMa_Truncated = SIGMA(1:truncated_dim,1:truncated_dim);
sigMa_Truncated = diag(sigMa_Truncated(sigMa_Truncated(:)~=0).^-0.5);
V = LX'*U(:,1:truncated_dim)*sigMa_Truncated;
W = LY'*Z(:,1:truncated_dim)*sigMa_Truncated;
disp(size(V))
disp(size(W))
q0   = W'*X0;
disp(q0)
% disp("V")
% disp(V)
% disp("W")
% disp(W)
%% simulate nonlinear reduced dynamics
nonode = @(t,p) auto_red_dyn(p,autData);
% option = odeset('RelTol',1e-8,'AbsTol',1e-10);
option = odeset('RelTol',1e-6,'AbsTol',1e-8);
% forward simulation of reduced dynamics
[~,pt] = ode45(nonode,tspan,p0,option);
% evaluation of W(p), bQ(t), and bM(t)
pt    = transpose(pt);
zauto = reduced_to_full_traj([],pt,Wmap);
bQ    = 2*ep*transpose(Q*V)*zauto;
bM    = 2*ep*transpose(Mhat*V)*zauto;
disp("bM")
disp(bM(:,end))
bt    = @(t) W'*Binv*Fext(t);


%% backward Riccati equation and compensation s(t)
% setup some matrices
Q2     = ep^2*transpose(V)*Q*V;
M2     = ep^2*transpose(V)*Mhat*V;
Lamhat = W'*Atil*V;
Bhat   = W'*Btil;
% backward simualtion of Riccati ODE
timeRt = tic;
% dPfun  = @(t,x) -RiccatiEqs(t,x,Lamhat,Bhat,Q2,Rhat); % minus sign for backward simulation
% [~,Pt] = ode15s(dPfun,tspan,M2(:),option);
% Pt     = flipud(Pt); % flip time order

dPfun  = @(t,x) RiccatiEqs(t,x,Lamhat,Bhat,Q2,Rhat); % minus sign for backward simulation
[~,Pt] = ode15s(dPfun,fliplr(tspan),M2(:),option);      % backward，P(t1) = M2, but in this case, we don't obtain what P(t1) is ? ? ?
Pt     = flipud(Pt);

fprintf('Time for backward simulation of Riccati ODE is %d\n',toc(timeRt));
% backward simulation of ODE for s
times  = tic;
% Ptfun  = @(t) interp1(tspan,Pt,t);
Ptfun  = griddedInterpolant(tspan,Pt); 
% bQt    = @(t) (interp1(tspan,transpose(bQ),t)).';
bQtmp  = griddedInterpolant(tspan,transpose(bQ)); bQt = @(t) bQtmp(t).';
Ptau   = @(tau) Ptfun(t1-tau); % transformation: t=t0+s=t1-tau
bQtau  = @(tau) bQt(t1-tau);
btau   = @(tau) bt(t1-tau);
dsfun  = @(tau,x) -compsEqs(tau,x,Ptau,Bhat,Rhat,Lamhat,bQtau,btau);
[~,st] = ode15s(dsfun,tspan-t0,-bM(:,end),option);
st     = flipud(st); % flip time order
% stfun  = @(t) (interp1(tspan,st,t)).';
stfunt = griddedInterpolant(tspan,st); stfun =  @(t) stfunt(t).';
fprintf('Time for backward simulation of compensated ODE is %d\n',toc(times));

%% forward simulation of qdot to yield q(t)
timeqt = tic;
odeq   = @(t,q) Lamhat*q+Bhat*closed_loop_control_lqr(t,q,Rhat,Bhat,Ptfun,stfun)+bt(t);
[~,qt] = ode15s(odeq,tspan,q0,option);  qt = qt.';
% evalute control input
ut     = closed_loop_control_lqr(tspan,qt,Rhat,Bhat,Ptfun,stfun); 
% evalute state of original system
zt     = zauto+ep*real(V*qt);
outdof = DS.Options.outDOF;
zout   = zt(outdof,:);
fprintf('Time for forward simulation of ODEs for modal coordinates is %d\n',toc(timeqt));
% invariance PDE error
res   = invariance_PDE_residual(DS,tspan,zt,ut);
% ratio = ratio_ext2int(DS,tspan,zt,ut);
%% output
traj      = struct();
traj.time = tspan;
traj.ut = ut; traj.pt = pt;
traj.qt = qt; traj.zt = zout; traj.res = res;

end