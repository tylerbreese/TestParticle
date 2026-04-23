A = readtable("sim_data_2026-04-22_13-51-33.csv");
sim = table2array(A);
B = readtable("split_data_2026-04-22_13-51-33.csv");
split = table2array(B);

En0 = sim(:,4) ./ 1.6e-9; 
En1 = sim(:,5) ./ 1.6e-9;
Ens = split(:,6) ./ 1.6e-9;
bins = logspace(-1,4,25);

n0 = histcounts(En0,bins);
n1 = histcounts(En1,bins);
ns = histcounts(Ens,bins);
n0(end+1) = 0.0; n1(end+1) = 0.0; ns(end+1) = 0.0;

figure()
hold on
grid on
scatter(bins,n0./sum(n0),'filled')
scatter(bins,n1./sum(n1),'filled')
scatter(bins,ns./sum(ns),'filled')
set(gca,'XScale','log','YScale','log');