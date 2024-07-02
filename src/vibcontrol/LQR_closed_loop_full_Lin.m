function traj = LQR_closed_loop_full_Lin(DS,z0,proj,tspan,autData,Wmap,masterModes,cont)
% LQR_closed_loop_full_Lin This routune presents LQR optimal control via SSM-based
% model reduction. Here we add closed-loop control strategy. For the
% perturbed linear part, we do not perform model reduction.

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

%% simulate nonlinear reduced dynamics
nonode = @(t,p) auto_red_dyn(p,autData);
option = odeset('RelTol',1e-8,'AbsTol',1e-10);
% forward simulation of reduced dynamics
[~,pt] = ode45(nonode,tspan,p0,option);
% evaluation of W(p), bQtil(t), and bMtil(t)
pt    = transpose(pt);
zauto = reduced_to_full_traj([],pt,Wmap);
bQtil  = 2*ep*transpose(Q)*zauto;
Binv   = DS.B\eye(DS.N);
bMtil  = 2*ep*transpose(Mhat)*zauto;
bttil  = @(t) Binv*Fext(t);

%% backward Riccati equation and compensation s(t)
% setup some matrices
Q2til  = ep^2*Q;
M2til  = ep^2*Mhat;
Atil   = Binv*DS.A;
Btil   = Binv*Bext;
% backward simualtion of Riccati ODE
timeRt = tic;
dPfun  = @(t,x) -RiccatiEqs(t,x,Atil,Btil,Q2til,Rhat); % minus sign for backward simulation
[~,Pt] = ode15s(dPfun,tspan,M2til(:),option);
Pt     = flipud(Pt); % flip time order
fprintf('Time for backward simulation of Riccati ODE is %d\n',toc(timeRt));
% backward simulation of ODE for s
times  = tic;
% Ptfun  = @(t) interp1(tspan,Pt,t);
Ptfun  = griddedInterpolant(tspan,Pt); 
% bQt    = @(t) (interp1(tspan,transpose(bQ),t)).';
bQtmp  = griddedInterpolant(tspan,transpose(bQtil)); bQt = @(t) bQtmp(t).';
Ptau   = @(tau) Ptfun(t1-tau); % transformation: t=t0+s=t1-tau
bQtau  = @(tau) bQt(t1-tau);
btau   = @(tau) bttil(t1-tau);
dsfun  = @(tau,x) -compsEqs(tau,x,Ptau,Btil,Rhat,Atil,bQtau,btau);
[~,st] = ode15s(dsfun,tspan-t0,-bMtil(:,end),option);
st     = flipud(st); % flip time order
% stfun  = @(t) (interp1(tspan,st,t)).';
stfunt = griddedInterpolant(tspan,st); stfun =  @(t) stfunt(t).';
fprintf('Time for backward simulation of compensated ODE is %d\n',toc(times));

%% forward simulation of X0dot to yield X0(t)
timeqt = tic;
odeq   = @(t,X) Atil*X+Btil*closed_loop_control_lqr(t,X,Rhat,Btil,Ptfun,stfun)+bttil(t);
[~,X0t] = ode15s(odeq,tspan,X0,option);  X0t = X0t.';
% evalute control input
ut     = closed_loop_control_lqr(tspan,X0t,Rhat,Btil,Ptfun,stfun); 
% evalute state of original system
zt     = zauto+ep*X0t;
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
traj.X0t = X0t; traj.zt = zout; traj.res = res;

end