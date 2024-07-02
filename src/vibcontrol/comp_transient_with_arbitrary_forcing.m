function traj = comp_transient_with_arbitrary_forcing(DS,z0,proj,tspan,autData,Wmap,masterModes,varargin)
% COMP_TRANSIENT_WITH_ARBITRARY_FORCING This function computes transient
% responses of high-dimensional systems via two-levels of reduction.

%% setup initial conditions
p0 = get_initial_p0(DS,masterModes,z0,proj,autData,Wmap);
X0 = (z0-reduced_to_full(p0,Wmap,[],0))/DS.epsilon;

%% simulate nonlinear reduced dynamics
nonode = @(t,p) auto_red_dyn(p,autData);
option = odeset('RelTol',1e-8,'AbsTol',1e-10);
outdof = DS.Options.outDOF;
% forward simulation of reduced dynamics
[~,pt] = ode45(nonode,tspan,p0,option);
% mapping it back to physical domain
pt   = transpose(pt);
z    = reduced_to_full_traj([],pt,Wmap);
znon = z(outdof,:);

%% simulate linear reduced dynamics
if numel(varargin)>0 && isnumeric(varargin{1})
    linModes = varargin{1};
    Uhat  = DS.spectrum.W(:,linModes);
    Vhat  = DS.spectrum.V(:,linModes);
    q0    = Uhat'*DS.B*X0;
    LamHat = diag(DS.spectrum.Lambda(linModes));
    linode = @(t,q) LamHat*q+Uhat'*DS.evaluate_Fext(t)/DS.epsilon;
    [~,qt] = ode45(linode,tspan,q0);
    qt     = transpose(qt);
    zlin   = Vhat*qt;
    zlin   = real(zlin(outdof,:));
else
    % call time integration transient to yield linear response 
    DSnonauto = DynamicalSystem();
    set(DSnonauto,'M',DS.M,'C',DS.C,'K',DS.K,'E',DS.E,'D',DS.D,'u',DS.u);
    om = 2*pi/tspan(end); nsteps = numel(tspan)-1;
    int_scheme = 'ode45';
    if numel(varargin)>0 && ischar(varargin{1}); int_scheme = varargin{1}; end
    [~, zlin] = time_integration_transient(DSnonauto,om,'nCycles',...
        1, 'nSteps', nsteps,...
        'integrationMethod',int_scheme,'outdof',outdof,'init',X0);
    zlin = zlin';
end

%% summation of nonlinear and linear parts
traj = struct();
traj.time = tspan;
traj.p = pt;
traj.z = znon+DS.epsilon*zlin;

end