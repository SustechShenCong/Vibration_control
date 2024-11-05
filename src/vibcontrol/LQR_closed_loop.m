function traj = LQR_closed_loop(DS,z0,proj,tspan,autData,Wmap,masterModes,linModes,cont)
% LQR_closed_loops This routune presents LQR optimal control via SSM-based
% model reduction. Here we add closed-loop control strategy.

%% setup of system and control
Q    = cont.Q;
Rhat = cont.Rhat;
Mhat = cont.Mhat;
ep   = DS.epsilon;
n    = size(DS.M,1);
Bext = [DS.D; zeros(size(DS.D))];
t0   = tspan(1); t1 = tspan(end);
if isempty(DS.E)
    Fext = @(t) zeros(2*n,1);
else
    Fext = @(t) [DS.E(t); zeros(n,1)];
end

%% setup initial conditions

p0   = get_initial_p0(DS,masterModes,z0,proj,autData,Wmap);
X0   = (z0-reduced_to_full(p0,Wmap,[],0))/ep;

Vhat = DS.spectrum.V(:,linModes);
Uhat = DS.spectrum.W(:,linModes);
q0   = Uhat'*DS.B*X0;


%% simulate nonlinear reduced dynamics
nonode = @(t,p) auto_red_dyn(p,autData);
option = odeset('RelTol',1e-8,'AbsTol',1e-10);
% forward simulation of reduced dynamics
% tic 
[~,pt] = ode45(nonode,tspan,p0,option);
% toc
% evaluation of W(p), bQ(t), and bM(t)
pt    = transpose(pt);
zauto = reduced_to_full_traj([],pt,Wmap);
bQ    = 2*ep*transpose(Q*Vhat)*(zauto-0);
bM    = 2*ep*transpose(Mhat*Vhat)*(zauto-0);
bt    = @(t) Uhat'*Fext(t);

%% backward Riccati equation and compensation s(t)
% setup some matrices
Q2     = ep^2*transpose(Vhat)*Q*Vhat;
M2     = ep^2*transpose(Vhat)*Mhat*Vhat;
Lamhat = diag(DS.spectrum.Lambda(linModes));
Bhat   = Uhat'*Bext;
% backward simualtion of Riccati ODE
timeRt = tic;
% dPfun  = @(t,x) -RiccatiEqs(t,x,Lamhat,Bhat,Q2,Rhat); % minus sign for backward simulation
% [~,Pt] = ode15s(dPfun,tspan,M2(:),option);
% Pt     = flipud(Pt); % flip time order

dPfun  = @(t,x) RiccatiEqs(t,x,Lamhat,Bhat,Q2,Rhat); % minus sign for backward simulation
% tic 
[~,Pt] = ode15s(dPfun,fliplr(tspan),M2(:),option);      % backward，P(t1) = M2, but in this case, we don't obtain what P(t1) is ? ? ?
% toc
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
% tic
[~,st] = ode15s(dsfun,tspan-t0,-bM(:,end),option);
% toc
st     = flipud(st); % flip time order
% stfun  = @(t) (interp1(tspan,st,t)).';
stfunt = griddedInterpolant(tspan,st); stfun =  @(t) stfunt(t).';
fprintf('Time for backward simulation of compensated ODE is %d\n',toc(times));

%% forward simulation of qdot to yield q(t)
timeqt = tic;
odeq   = @(t,q) Lamhat*q+Bhat*closed_loop_control_lqr(t,q,Rhat,Bhat,Ptfun,stfun)+bt(t);
tic 
[~,qt] = ode15s(odeq,tspan,q0,option);  qt = qt.';
% evalute control input
ut     = closed_loop_control_lqr(tspan,qt,Rhat,Bhat,Ptfun,stfun); 
toc
% evalute state of original system
zt     = zauto+ep*real(Vhat*qt);
outdof = DS.Options.outDOF;
zout   = zt(outdof,:);
zout_auto   = zauto(outdof,:);
fprintf('Time for forward simulation of ODEs for modal coordinates is %d\n',toc(timeqt));
% invariance PDE error
% res   = invariance_PDE_residual(DS,tspan,zt,ut);
% ratio = ratio_ext2int(DS,tspan,zt,ut);
%% output
traj      = struct();
traj.time = tspan;
traj.ut = ut; traj.pt = pt;
traj.qt = qt; traj.zt = zout; traj.zt_auto = zout_auto;
traj.res = [];

end