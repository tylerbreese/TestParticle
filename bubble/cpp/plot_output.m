A = readtable("sim_data_2026-04-22_20-48-44");
sim = table2array(A);
B = readtable("split_data_2026-04-22_20-48-44.csv");
split = table2array(B);

x0 = sim(:,2) ./ 1e5;
xf = sim(:,3) ./ 1e5;
En0 = sim(:,4) ./ 1.6e-9; 
En1 = sim(:,5) ./ 1.6e-9;
Ens = split(:,6) ./ 1.6e-9;
bins = logspace(-1,4,25);
%%
n0 = histcounts(En0,bins);
n1 = histcounts(En1,bins);
ns = histcounts(Ens,bins);
n0(end+1) = 0.0; n1(end+1) = 0.0; ns(end+1) = 0.0;
%%
figure()
hold on
grid on
scatter(bins,n0./sum(n0),'filled')
scatter(bins,n1./sum(n1),'filled')
%scatter(bins,ns./sum(ns),'filled')
set(gca,'XScale','log','YScale','log');

figure()
hold on
grid on
scatter(x0,En0,'filled')
scatter(xf,En1,'filled')
set(gca,'XScale','log','YScale','log');