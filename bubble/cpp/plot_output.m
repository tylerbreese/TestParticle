A = readtable("output/sim_data_2026-05-03_12-13-07.csv");
sim = table2array(A);
B = readtable("output/sim_data_2026-05-03_12-13-07.csv");
split = table2array(B);

x0 = sim(:,2) ./ 1e5;
xf = sim(:,3) ./ 1e5;
xs = split(:,4) ./ 1e5;
En0 = sim(:,4); 
En1 = sim(:,5);
Ens = split(:,6);

%%
% bins = logspace(0,1,10);
% n0 = histcounts(En0,bins);
% n1 = histcounts(En1,bins);
% ns = histcounts(Ens,bins);
% n0(end+1) = 0.0; n1(end+1) = 0.0; ns(end+1) = 0.0;
% 
% figure()
% hold on
% grid on
% scatter(bins,n0./sum(n0),'filled')
% scatter(bins,n1./sum(n1),'filled')
% scatter(bins,ns./sum(ns),'filled')
% set(gca,'XScale','log','YScale','log');

[n0,bins0] = histcounts(En0,25);
[n1,bins1] = histcounts(Enf,25);
n0(end+1) = 0.0; n1(end+1) = 0.0;
figure()
hold on
grid on
scatter(bins0,n0./sum(n0))
scatter(bins1,n1./sum(n1))
set(gca,'XScale','log','YScale','log');

figure()
hold on
grid on
scatter(x0,En0,'filled')
scatter(xf,En1,'filled')
scatter(xs,Ens,'filled')
set(gca,'YScale','log');