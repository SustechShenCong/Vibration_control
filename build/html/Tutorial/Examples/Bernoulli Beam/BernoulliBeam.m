%% Euler Bernoulli beam with cubic spring and damper
% In this notebook, we consider the forced response curve of a Euler Bernoulli 
% beam with cubic spring and damper at the free end, subject to a harmonic base 
% excitation. In the frame attached to the oscillating base, the beam is subject 
% to harmonic excitation in the form $m\Omega^2\cos\Omega t$, where $\Omega$ is 
% the excitation frequency.
%% Generate model

clear all
nElements = 5; % number of beam elements
kappa = 4; % cubic spring coefficient
gamma = 0; % cubic damping coefficient
[M,C,K,fnl,~] = build_model(kappa, gamma, nElements);
n = length(M);
%% Dynamical system setup 
% We consider the forced system
% 
% $$\mathbf{M}\ddot{\mathbf{x}}+\mathbf{C}\dot{\mathbf{x}}+\mathbf{K}\mathbf{x}+\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})=\epsilon\mathbf{M}\cdot\mathbf{1}(\Omega^2\cos{\Omega}t),$$
%
% In our code, the dynamical system is represented by the system matrices
% and an implementation of the nonlinear internal forces.
%

DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
set(DS.Options,'BaseExcitation',true); % Frequency dependent forcing amplitude
%% 
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{M}\cdot\mathbf{1}(\Omega^2\cos{\Omega}t)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
%
% where $\mathbf{f}_0 =  \mathbf{M}\cdot\mathbf{1}$.
%

epsilon = 5e-5;
f_0 = M*ones(n,1);
%% 
%
% The Fourier coefficients of the forcing and the forcing amplitude need to
% be passed to the dynamical system.
%

kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
%% Linear Modal analysis and SSM setup
%
% Now the reduced spectral content linear part of the dynamicals system is
% analysed. For this a set number (see DSOptions) of modes is computed. 
%

[V,D,W] = DS.linear_spectral_analysis();
%% 
%
% *Choose Master subspace (perform resonance analysis)*
%
% Next the spectral subspace over which the invariant SSM is to be
% constructed is chosen. The notation in which the SSM is to be computed is
% set to the multi-index notation. Relative tolerances are triggered at a
% relative modal distance of 0.1 (see DSOptions), which allows us to get an overview of how separated the modes are prior to computing the SSM.
% The spectral subspace consisting of the primary mode-pair is chosen as a base for the invariant
% manifold.
%


S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2]; 
S.choose_E(masterModes);
%% Forced response curves using SSMs
%
% Obtaining *forced response curve* in reduced-polar coordinate
%
% Now we setup the SSM computation, to compute the forced response of the
% beam, if driven in resonance with the first mode.
%
% The SSM is to be computed up to order 3 and 5, to monitor the convergence
% of the parametrisation and reduced dynamics expansions. We wish to
% display the resulting FRCs along the output DOFs at the tip of the beam
% and set the corresponding parameters accordingly.

order = [3,5];      % SSM approximation order
outdof = [n-2 n-1]; % degree of freedom at which output is displayed
%% 
% Now we setup the options for the computation of the SSM. For the
% computation of the FRC we choose the level set method. The steady states
% are therefore computed directly from the polar ODE of the reduced
% dynamics (see <FRC computation>).
%

set(S.FRCOptions, 'method','level set')  

%%
%
% We set the relative tolerance for internal resonances and choose to include nonautonomous 
% contributions to the invariant manifold and reduced dynamics. For the
% mesh over which the polar ODE is evaluated to find fixedpoints we choose
% a 200 by 100 grid in the polar coordinates. The number of frequencies on which the ODE should be 
% evaluated is set to 100. The parameter rhoScale is used to scale the
% amplituderange in which the reduced system response is to be analysed.
%

set(S.Options, 'reltol', 1,'IRtol',0.02,'contribNonAuto',true)
set(S.FRCOptions, 'nRho', 200,  'nPsi', 100, 'nPar', 100,'rhoScale', 2 )

%% 
% We now choose the frequency range over which the FRC is to be computed.
%

omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.98 1.05];
%% 
% Finally extract forced response curve via a computation of the SSM and
% the reduced dynamics on it.

FRC = S.extract_FRC('freq',omegaRange,order);