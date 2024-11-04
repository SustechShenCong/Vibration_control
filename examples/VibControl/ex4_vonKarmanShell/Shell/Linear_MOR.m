function res = Linear_MOR(DS,linModes)
% LQR_closed_loops This routune presents LQR optimal control via SSM-based
% model reduction. Here we add closed-loop control strategy.

%% By Hankle singular value
Vhat = DS.spectrum.V(:,1:2*linModes);
Uhat = DS.spectrum.W(:,1:2*linModes);
Lamhat = DS.spectrum.Lambda(1:2*linModes);

outdof = DS.Options.outDOF;
m      = length(outdof);
n      = size(DS.M,1);
HSV     = zeros(linModes,1);
DCGain_ = zeros(linModes,1);
Bext = [DS.D; zeros(size(DS.D))];

for k = 1:linModes
    W = Uhat(:, 2*k-1:2*k);
    V = Vhat(:, 2*k-1:2*k);
    A = diag(Lamhat(2*k-1:2*k));
    B = W' * Bext;
    % C        = eye(size(A)); % two kind of observation matrix C
    C        = zeros(m, n);
    for i = 1:m
        C(i,outdof(i)) = 1;
    end
    C        = [C,   zeros(size(C))];
    C        = C * V;
    D        = zeros(length(outdof),size(B,2));
    sys      = ss(A,B,C,D);
    
    LX = lyap(A,B*B');
    LY = lyap(A',C'*C);

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

HSV_per = zeros(size(HSV));
HSV_sum = sum(HSV);
HSV_per(1) = HSV(1) / HSV_sum;
for i = 2:length(HSV)
    % HSV_per(i) = HSV_per(i-1) + HSV(i) / HSV_sum;
    HSV_per(i) = HSV(i) / HSV_sum;
end


DC_per = zeros(size(DCGain_));
DC_sum = sum(DCGain_);
DC_per(1) = DCGain_(1) / DC_sum;
for i = 2:length(DCGain_)
    % DC_per(i) = DC_per(i-1) + DCGain_(i) / DC_sum;
    DC_per(i) = DCGain_(i) / DC_sum;
end

figure; hold on
subplot(1,2,1)
% plot(HSV_per,'k-',"LineWidth",1)
bar(HSV_per,'Linewidth',0.5,"FaceColor","b",'Edgecolor',"none","FaceAlpha",0.85)
% plot([0, length(HSV_per)], [0.95, 0.95], 'k--')
xlabel('Number of mode pair','Interpreter','latex')
zk = strcat('Normalized MHSV');
% zk = strcat('$Displacement Propotion$');
ylabel(zk,'Interpreter','latex');
% title('HSV of each DOFs - 2D reduction')
set(gca,'FontSize',19.8);
grid on, axis tight
% legend('Hankle Singular Value','Interpreter',"latex"); legend boxoff
box on
ylim([0,1])
xticks([1, 4, 7, 10])
set(gca, 'LineWidth', 2);
% set(gcf, 'Position', [0 0 600 500]);
% print('-depsc', 'OscillatorChain_ControlPolicy_outdof5_woFeedback.eps')
% print -djpeg -r300 OscillatoChain_HSVs.jpg;



% figure; hold on
subplot(1,2,2)
% plot(DC_per,'k-',"LineWidth",1)
bar(DC_per,'Linewidth',0.5,"FaceColor",	"b",'Edgecolor',"none","FaceAlpha",0.85)
% plot([0, length(HSV_per)], [0.95, 0.95], 'k--')
xlabel('Number of mode pair','Interpreter',"latex")
zk = strcat('Normalized DCgain');
% zk = strcat('$Displacement Propotion$');
ylabel(zk,'Interpreter','latex');
% title('DCGain of each DOFs - 2D reduction')
set(gca,'FontSize',19.8);
grid on, axis tight
% legend('DCGain','Interpreter',"latex"); legend boxoff
ylim([0,1])
xticks([1, 4, 7, 10])
box on

set(gca, 'LineWidth', 2);
set(gcf, 'Position', [0 0 700 500]);
% print('-depsc', 'OscillatorChain_ControlPolicy_outdof5_woFeedback.eps')
% print -djpeg -r300 OscillatoChain_DCGains.jpg;


disp('MHSV sum of the first 2 pairs')
disp(sum(HSV_per(1:2)))
disp('DCgain sum of the first 2 pairs')
disp(sum(DC_per(1:2)))
res{1} = HSV_per;
res{2} = DC_per;


end
