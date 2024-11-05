function res = Linear_MOR(DS,linModes, phiend)
% LQR_closed_loops This routune presents LQR optimal control via SSM-based
% model reduction. Here we add closed-loop control strategy.

%% By Hankle singular value
n = size(DS.M,1);
Vhat = DS.spectrum.V(:,1:2*linModes);
Uhat = DS.spectrum.W(:,1:2*linModes);
Lamhat = DS.spectrum.Lambda(1:2*linModes);

HSV     = zeros(linModes,1);
DCGain_ = zeros(linModes,1);
Bext = [DS.D; zeros(size(DS.D))];

for k = 1:linModes
    W = Uhat(:, 2*k-1:2*k);
    V = Vhat(:, 2*k-1:2*k);
    A = diag(Lamhat(2*k-1:2*k));
    B = W' * Bext;
    C = [phiend; zeros(n,1)]';
    C        = C * V;

    D        = 0;
    sys      = ss(A,B,C,D);
    LX = lyap(A,  B*B');
    LY = lyap(A', C'*C);

    [U, SIGMA, Z] = svd(LX*LY');
    SIGMA = sqrt(SIGMA);

    HSV(k)     = max(diag(SIGMA));
    gains      = dcgain(sys);
    DCGain_(k) = norm(gains,2);
end

HSV(1) = [];
HSV_per = zeros(size(HSV));
HSV_sum = sum(HSV);
HSV_per(1) = HSV(1) / HSV_sum;
for i = 2:length(HSV)
    HSV_per(i) = HSV(i) / HSV_sum;
end

DCGain_(1) = [];
DC_per = zeros(size(DCGain_));
DC_sum = sum(DCGain_);
for i = 1:length(DCGain_)
    DC_per(i) = DCGain_(i) / DC_sum;
end

figure; hold on
subplot(1,2,1)
bar([2,3,4], HSV_per,'Linewidth',0.5,"FaceColor","b",'Edgecolor',"none","FaceAlpha",0.85)
xlabel('Number of mode pair','Interpreter','latex')
zk = strcat('Normalized MHSV');
ylabel(zk,'Interpreter','latex');
grid on, axis tight
set(gca,'FontSize',18); grid on, axis tight
box on
ylim([0,1])
xticks([2, 3, 4])
set(gca, 'LineWidth', 2);


subplot(1,2,2)
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

res{1} = HSV_per;
res{2} = DC_per;

end
