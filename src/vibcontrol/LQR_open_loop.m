function traj = LQR_open_loop(DS,z0,proj,tspan,autData,Wmap,masterModes,linModes,cont)
% LQR_open_loops This routune presents LQR optimal control via SSM-based
% model reduction

%% setup of system and control
Q    = cont.Q;
Rhat = cont.Rhat;
Mhat = cont.Mhat;
ep   = DS.epsilon;
n    = size(DS.M,1);
m    = numel(linModes);
Bext = [DS.D; zeros(size(DS.D))];
t0   = tspan(1);          t1 = tspan(end);
dt   = tspan(2)-tspan(1); T  = t1-t0;
if isempty(DS.E)
    Fext = @(t) zeros(2*n,1);
else
    Fext = @(t) [DS.E(t); zeros(n,1)];
end

%% setup initial conditions
p0 = get_initial_p0(DS,masterModes,z0,proj,autData,Wmap);
X0 = (z0-reduced_to_full(p0,Wmap,[],0))/ep;
Vhat  = DS.spectrum.V(:,linModes);
Uhat  = DS.spectrum.W(:,linModes);
disp("V")
disp(Vhat)
disp(size(Vhat))
q0    = Uhat'*DS.B*X0;

%% simulate nonlinear reduced dynamics
nonode = @(t,p) auto_red_dyn(p,autData);
option = odeset('RelTol',1e-8,'AbsTol',1e-10);
% forward simulation of reduced dynamics
[~,pt] = ode45(nonode,tspan,p0,option);
% evaluation of W(p), bQ(t), and bM(t)
pt    = transpose(pt);
zauto = reduced_to_full_traj([],pt,Wmap);
bQ    = 2*ep*transpose(Q*Vhat)*zauto;
bM    = 2*ep*transpose(Mhat*Vhat)*zauto;
bt    = @(t) Uhat'*Fext(t);

%% compute fundamental matrix and Duhamel term
% setup some matrices
Q2 = ep^2*transpose(Vhat)*Q*Vhat;
M2 = ep^2*transpose(Vhat)*Mhat*Vhat;
Lamhat = diag(DS.spectrum.Lambda(linModes));
Bhat = Uhat'*Bext;
Hmat = [Lamhat 0.5*Bhat*(Rhat\transpose(Bhat)); 2*Q2 -Lamhat];
% transition matrix from t0 to t1
phiT  = expm(Hmat*T); % or precise integration via exp_mat(Hmat,T)
% Duhamel term
bfun   = @(t) [bt(t); (interp1(tspan,transpose(bQ),t)).'];
odefun = @(t,x) Hmat*x+bfun(t);
[~,Duham] = ode45(odefun,tspan,zeros(2*m,1),option); Duham = Duham.';
% or precise integration via 
% Duham = pim_linf(Hmat,bfun,zeros(2*m,1),dt,t0,t1);
% compute initial costate
d1t   = Duham(1:m,:);      d2t   = Duham(m+1:2*m,:);
phi11 = phiT(1:m,1:m);     phi12 = phiT(1:m,m+1:2*m);
phi21 = phiT(m+1:2*m,1:m); phi22 = phiT(m+1:2*m,m+1:2*m);
mub0  = -(phi22+2*M2*phi12)\(bM(:,end)+phi21*q0+2*M2*(phi11*q0+d1t(:,end))+d2t(:,end));

%% evaluate control input and state
% forward simulation with obtained initial condition
[~,qmubt] = ode45(odefun,tspan,[q0;mub0],option); qmubt = qmubt.';
% or precise integration via 
% qmubt = pim_linf(Hmat,bfun,[q0;mub0],dt,t0,t1);
qt    = qmubt(1:m,:);
mubt  = qmubt(m+1:2*m,:);

% evalute control input
ut    = real(0.5*(Rhat\transpose(Bhat)*mubt)); 

% evalute state of original system
zt     = zauto+ep*real(Vhat*qt);
outdof = DS.Options.outDOF;
zout   = zt(outdof,:);

%% output
traj      = struct();
traj.time = tspan;
traj.ut = ut; traj.pt = pt;
traj.qt = qt; traj.zt = zout;

end