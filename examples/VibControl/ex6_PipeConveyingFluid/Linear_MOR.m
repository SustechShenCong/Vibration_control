function res = Linear_MOR(DS,linModes, phiend)
% LQR_closed_loops This routune presents LQR optimal control via SSM-based
% model reduction. Here we add closed-loop control strategy.

%% By Hankle singular value
n = size(DS.M,1);
Vhat = DS.spectrum.V(:,1:2*linModes);
Uhat = DS.spectrum.W(:,1:2*linModes);
Lamhat = DS.spectrum.Lambda(1:2*linModes);

% outdof = DS.Options.outDOF;

HSV     = zeros(linModes,1);
DCGain_ = zeros(linModes,1);
Bext = [DS.D; zeros(size(DS.D))];

for k = 1:linModes
    W = Uhat(:, 2*k-1:2*k);
    V = Vhat(:, 2*k-1:2*k);
    A = diag(Lamhat(2*k-1:2*k));
    B = W' * Bext;
    % C        = eye(size(A)); % two kind of observation matrix C
    % C        = zeros(length(outdof),size(DS.M,1));
    % C(1,outdof(1)) = 1;
    % C(2,outdof(2)) = 1;
    % for Pipe case, the end displacement is our interest, and we need to 
    % C        = [C,   zeros(size(C))];
    C = [phiend; zeros(n,1)]';
    C        = C * V;

    D        = 0;
    sys      = ss(A,B,C,D);
    LX = lyap(A,  B*B');
    LY = lyap(A', C'*C);

    % LX = lyapchol(A,B);
    % LY = lyapchol(A',C');
    [U, SIGMA, Z] = svd(LX*LY');
    SIGMA = sqrt(SIGMA);

    HSV(k)     = max(diag(SIGMA));
    gains      = dcgain(sys);
    % DCGain_(k) = abs(gains(1)); % have 2 gians, displacement and velocity, keep only displacement
    % DCGain_(k) = gains(1); % have 2 gians, displacement and velocity
    DCGain_(k) = norm(gains,2);
end

HSV(1) = [];
HSV_per = zeros(size(HSV));
HSV_sum = sum(HSV);
HSV_per(1) = HSV(1) / HSV_sum;
for i = 2:length(HSV)
    % HSV_per(i) = HSV_per(i-1) + HSV(i) / HSV_sum;
    HSV_per(i) = HSV(i) / HSV_sum;
end

DCGain_(1) = [];
DC_per = zeros(size(DCGain_));
DC_sum = sum(DCGain_);
% DC_per(1) = DCGain_(1) / DC_sum;
for i = 1:length(DCGain_)
    % DC_per(i) = DC_per(i-1) + DCGain_(i) / DC_sum;
    DC_per(i) = DCGain_(i) / DC_sum;
end

figure; hold on
subplot(1,2,1)
% plot(HSV_per,'k-',"LineWidth",1)
bar([2,3,4], HSV_per,'Linewidth',0.5,"FaceColor","b",'Edgecolor',"none","FaceAlpha",0.85)
% plot([0, length(HSV_per)], [0.95, 0.95], 'k--')
xlabel('Number of mode pair','Interpreter','latex')
zk = strcat('Normalized MHSV');
% zk = strcat('$Displacement Propotion$');
ylabel(zk,'Interpreter','latex');
% title('HSV of each DOFs - 2D reduction')
grid on, axis tight
% legend('Hankle Singular Value','Interpreter',"latex"); legend boxoff
set(gca,'FontSize',18); grid on, axis tight
box on
ylim([0,1])
xticks([2, 3, 4])
set(gca, 'LineWidth', 2);
% set(gcf, 'Position', [0 0 500 400]);
% % print('-depsc', 'OscillatorChain_outdof5_woFeedback.eps')
% print -djpeg -r300 Fluid_MHSVs2.jpg;


% figure; hold on
subplot(1,2,2)
% plot(DC_per,'k-',"LineWidth",1)
bar([2,3,4], DC_per,'Linewidth',0.5,"FaceColor",	"b",'Edgecolor',"none","FaceAlpha",0.85)
xlabel('Number of mode pair','Interpreter',"latex")
zk = strcat('Normalized DCgain');
ylabel(zk,'Interpreter','latex');
grid on, axis tight
set(gca,'FontSize',18); grid on, axis tight
ylim([0,1])
box on
xticks([2, 3, 4])
set(gca, 'LineWidth', 2);
set(gcf, 'Position', [0 0 700 450]);
% print('-depsc', 'OscillatorChain_outdof5_woFeedback.eps')
% print -djpeg -r300 Fluid_DCGains2.jpg;

res{1} = HSV_per;
% res{1} = [];
res{2} = DC_per;
% res = SIGMA(1:2:end, 1:2:end);
% X0(1321:end) = 0;


end
