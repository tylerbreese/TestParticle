%% read in ddata
h5file(1) = "output_AGU/run_09032025_VSW_400_r_2p5_s0p5.h5";
h5file(2) = "output_AGU/run_09032025_VSW_400_r_4_s0p5.h5";
% h5file(1) = "proton_2p5_s0p5.h5";
% h5file(2) = "proton_4_s0p5.h5";

X_2p5     = h5read(h5file(1),'/position');
V_2p5     = h5read(h5file(1),'/velocity');
Split_2p5 = h5read(h5file(1),'/split');
Vmag_2p5  = vecnorm(V_2p5,2,2);

X_4     = h5read(h5file(2),'/position');
V_4     = h5read(h5file(2),'/velocity');
Split_4 = h5read(h5file(2),'/split');
Vmag_4  = vecnorm(V_4,2,2);

V_SW        = 400e5 .* [1.0; 0.0; 0.0;];
B0          = 0.05e-5 .* [0.0; 0.0; 1.0;];

U1 = vecnorm(V_SW,2,1);
B1 = vecnorm(B0,2,1);

particleData = h5read(h5file(1),'/particle');
m = particleData(1);
q = particleData(2);
Z = particleData(3);
A = particleData(4);

mask2p5 = X_2p5(:,1,2) > 0;
mask4   = X_4(:,1,2) > 0;

En_2p5 = 0.5 .* m .* squeeze(Vmag_2p5).^2 ./ 1.6e-6;
En_4   = 0.5 .* m .* squeeze(Vmag_4).^2 ./ 1.6e-6;

En_split_2p5 = 0.5 .* m .* vecnorm(Split_2p5,2,2).^2 ./ 1.6e-6;
En_split_4   = 0.5 .* m .* vecnorm(Split_4,2,2).^2 ./ 1.6e-6;

bins = logspace(-3,3,100);

n0 = histcounts(En_2p5(:,3), bins);
n_2p5 = histcounts([En_2p5(mask2p5,2); 0.5.*En_split_2p5], bins);
n_4   = histcounts([En_4(mask4,2); 0.5.*En_split_4], bins);

n_2p5(end+1) = 0.0;
n_4(end+1) = 0.0;
n0(end+1) = 0.0;

% I0   = trapz(bins,n0);
% I2p5 = trapz(bins,n_2p5);
% I4   = trapz(bins,n_4);
I0   = sum(n0);
I2p5 = sum(n_2p5);
I4   = sum(n_4);

%%
jE_0 = (n0) ./( (4*pi) .* (sqrt(bins)./sqrt(2/(A*938))) );
jE_2p5 = (n_2p5) ./( (4*pi) .* (sqrt(bins)./sqrt(2/(A*938))) );
jE_4 = (n_4) ./( (4*pi) .* (sqrt(bins)./sqrt(2/(A*938))) );

% dJ/dE = p^2 f
conv = 1.87e16; %g cm /s to MeV /c
p_2p5 = conv .* m .* squeeze(Vmag_2p5);
p_2p5_split = conv .* m .* vecnorm(Split_2p5,2,2);
p_4   = conv .* m .* squeeze(Vmag_4);
p_4_split = conv .* m .* vecnorm(Split_4,2,2);
c = 3e10;

np_0   = histcounts(p_2p5(:,3).^2,(2.*(A*938).*bins));
np_2p5 = histcounts([p_2p5(:,2).^2; p_2p5_split.^2],(2.*(A*938).*bins));
np_4   = histcounts([p_4(:,2).^2; p_4_split.^2],(2.*(A*938).*bins));
% np_2p5 = histcounts(p_2p5(:,2).^2,(2.*3752.*bins));
% np_4   = histcounts(p_4(:,2).^2,(2.*3752.*bins));

np_0(end+1) = 0.0;
np_2p5(end+1) = 0.0;
np_4(end+1) = 0.0;


%%
figure()
hold on
grid on
plot(bins,n_2p5./I2p5,'--b','LineWidth',2)
plot(bins,n_4./I4,'--r','LineWidth',2)
plot(bins,n0./I0,'-k','LineWidth',2)
xlabel('Energy (MeV)');
ylabel('Normalized Units');
legend('B_2 / B_1 = 2.5','B_2 / B_1 = 4.0','Initial');
title('Energy Distribution of Particles');
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
saveas(gcf,'energy','fig')
saveas(gcf,'energy','png')

%% flux
figure()
hold on
grid on
scatter(bins,jE_2p5,'b','filled')
scatter(bins,jE_4,'r','filled')
scatter(bins,jE_0,'k','filled')
ax=gca;
ax.YScale = 'log';
ax.XScale = 'log';
xlabel('Energy (MeV)');
ylabel('Flux (#/cm^2/s/sr/MeV)');
title('Particle Flux Distribution');
legend('B_2 / B_1 = 2.5','B_2 / B_1 = 4.0','Initial');
saveas(gcf,'flux','fig')
saveas(gcf,'flux','png')
%%
nH = 0.115;
%xs = 5e-15./sqrt(1e7.*bins);

a = [4.15, 0.531, 67.3];
xs = @(E) 1e-16 .* (a(1) - a(2) .* log(E)).^2 .* (1 - exp(-a(3)./(E))).^(4.5);
lambda = 1 ./ (xs(1000.*bins).*nH);
DR = 35 .* 1.496e13;

ENA_2p5 = (jE_2p5) .* (DR ./ lambda);
ENA_4   = (jE_4) .* (DR ./ lambda);

LECP = [4.558E+02  4.481E+00  2.147E+02  1.958E+00  1.030E+02  1.055E+00  4.924E+01  6.262E-01  1.771E+01  1.835E-01  1.261E+01  1.317E-01  5.923E+00  5.645E-02  1.358E+00  2.486E-02];
V2bins = [3.25*0.01,6*0.01, 0.1, 1.75*0.1, 3.25*0.1, 7*0.1, 1.5, 3];
%V2lambda = 1 ./ (nH .* 5e-15./sqrt(1e7.*V2bins));
V2lambda = 1 ./ (nH .* xs(1000.*V2bins));
V2ENA = LECP(1:2:16).* (DR ./ V2lambda);
%note energy bins might by slightly off
%check with joe 
% some colors
dg = hex2rgb('#006400');
pp = hex2rgb('#6F00BA');
or = hex2rgb('#FF9800');
bb = hex2rgb('#4472C4');

IBEX = [
    0.015, 1.15e8
    0.029, 3.5e7
    0.055, 1.35e7
    0.110, 1.05e6
    0.209, 4000
    0.439, 1000
    0.872, 380
    1.821, 125
    0.710,  500
    1.110, 300
    1.740, 121
    2.730, 60
    4.290, 14
    5+(13-5)/2,	12
    13+(24-13)/2, 0.5
    24+(35-24)/2, 0.06
    35+(55-35)/2, 0.02
    ];

figure()
hold on
grid on
scatter(bins.*1000,ENA_2p5,[],pp,'filled')
scatter(bins.*1000,ENA_4,[],or,'filled')
scatter(V2bins.*1000,V2ENA,[],dg,'square','filled')
scatter(IBEX(:,1),IBEX(:,2),[],bb,'square','filled')
ax=gca;
ax.YScale = 'log';
ax.XScale = 'log';
xlabel('Energy (keV)');
ylabel('Flux (#/cm^2/s/sr/keV)');
title('ENA Flux Distribution');
lgd = legend('B_2/B_1=2.5','B_2/B_1=4.0','H+ V2 LECP','IBEX+CASSINI', 'Location', 'northeast');

saveas(gcf,'ENA','fig')
saveas(gcf,'ENA','png')

%%
nuc = A;

figure()
hold on
grid on
scatter(bins./nuc,jE_2p5,'b','filled')
scatter(bins./nuc,jE_4,'r','filled')
scatter(bins./nuc,jE_0,'k','filled')
scatter(V2bins,LECP(1:2:16),[],dg,'square','filled')
ax=gca;
ax.YScale = 'log';
ax.XScale = 'log';
xlabel('Energy (MeV/nucleon)');
ylabel('Flux (#/cm^2/s/sr/MeV)');
title('Particle Flux Distribution');
legend('B_2 / B_1 = 2.5','B_2 / B_1 = 4.0','Initial','V2 LECP');
saveas(gcf,'flux1','fig')
saveas(gcf,'flux1','png')


