function res = Linear_MOR(DS,linModes)
% LQR_closed_loops This routune presents LQR optimal control via SSM-based
% model reduction. Here we add closed-loop control strategy.

%% By Hankle singular value
Vhat = DS.spectrum.V(:,1:2*linModes);
Uhat = DS.spectrum.W(:,1:2*linModes);
Lamhat = DS.spectrum.Lambda(1:2*linModes);

outdof = DS.Options.outDOF;

HSV     = zeros(linModes,1);
DCGain_ = zeros(linModes,1);
Bext = [DS.D; zeros(size(DS.D))];

for k = 1:linModes
    W = Uhat(:, 2*k-1:2*k);
    V = Vhat(:, 2*k-1:2*k);
    A = diag(Lamhat(2*k-1:2*k));
    B = W' * Bext;
    C        = zeros(length(outdof),size(DS.M,1));
    C(1,outdof(1)) = 1;
    C(2,outdof(2)) = 1;
    C        = [C,   zeros(size(C))];
    C        = C * V;
    D        = zeros(length(outdof),size(B,2));
    sys      = ss(A,B,C,D);
    
    LX = lyap(A,B*B');
    LY = lyap(A',C'*C);

    [U, SIGMA, Z] = svd(LX*LY');
    SIGMA = sqrt(SIGMA);
    
    HSV(k)     = max(diag(SIGMA));
    gains      = dcgain(sys);
    DCGain_(k) = norm(gains,2);
end

HSV_per = zeros(size(HSV));
HSV_sum = sum(HSV);
HSV_per(1) = HSV(1) / HSV_sum;
for i = 2:length(HSV)
    HSV_per(i) = HSV(i) / HSV_sum;
end


DC_per = zeros(size(DCGain_));
DC_sum = sum(DCGain_);
DC_per(1) = DCGain_(1) / DC_sum;
for i = 2:length(DCGain_)
    DC_per(i) = DCGain_(i) / DC_sum;
end

figure; hold on

bar(HSV_per,'Linewidth',0.5,"FaceColor","b",'Edgecolor',"none","FaceAlpha",0.85)
xlabel('Number of mode pair','Interpreter','latex')
zk = strcat('Normalized MHSV');
ylabel(zk,'Interpreter','latex');
set(gca,'FontSize',24);
grid on, axis tight
set(gca,'FontSize',24); grid on, axis tight
box on
ylim([0,1])

set(gca, 'LineWidth', 2);
set(gcf, 'Position', [0 0 600 500]);

figure; hold on
bar(DC_per,'Linewidth',0.5,"FaceColor",	"b",'Edgecolor',"none","FaceAlpha",0.85)
xlabel('Number of mode pair','Interpreter',"latex")
zk = strcat('Normalized DCgain');
ylabel(zk,'Interpreter','latex');
set(gca,'FontSize',24);
grid on, axis tight
set(gca,'FontSize',24); grid on, axis tight
ylim([0,1])
box on

set(gca, 'LineWidth', 2);
set(gcf, 'Position', [0 0 600 500]);


disp('MHSV sum of the first 2 pairs')
disp(sum(HSV_per(1:2)))
disp('DCgain sum of the first 2 pairs')
disp(sum(DC_per(1:2)))
res{1} = HSV_per;
res{2} = DC_per;


end