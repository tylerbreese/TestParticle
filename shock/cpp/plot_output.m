folder = "IPShock_90";
data = dir(folder);
data = {data.name};
data = data(3:end);

sim_data = data(contains(data,'sim_data'));
split_data = data(contains(data,'split_data'));

for ii = 1:length(sim_data)

    A = readtable(fullfile(folder,sim_data{ii}));
    sim = table2array(A);
    En0(:,ii) = sim(:,4) ./ 1.6e-9; 
    
    En1(:,ii) = sim(:,5) ./ 1.6e-9;
    % B = readtable(fullfile(folder,split_data{ii}));
    % split = table2array(B);
    %Ens(:,ii) = split(:,6) ./ 1.6e-9;
    

end
%%
En0 = reshape(En0,[],1);
En1 = reshape(En1,[],1);
bins = logspace(-1,6,100);

[n0,b0] = histcounts(En0,bins);
[n1,b1] = histcounts(En1,bins);
n0(end+1) = 0.0; n1(end+1) = 0.0; 
%Ens = reshape(Ens,[],1);
%ns = histcounts(Ens,bins);
%ns(end+1) = 0.0;


fig1 = figure();
hold on
grid on
% scatter(bins,n0,'filled')
% scatter(bins,n1,'filled')
scatter(bins,n0./sum(n0),'filled')
scatter(bins,n1./sum(n1),'filled')
%scatter(bins,ns./sum(ns),'filled')
set(gca,'XScale','log','YScale','log');
xlabel('Energy (keV)')
ylabel('Normalized Units')
title('He+ PUI Distribution Function')
saveas(fig1,'graph1.png')
%%
A = readtable("sim_data_2026-05-31_13-04-34.csv");
sim = table2array(A);
B = readtable("split_data_2026-05-31_13-04-34.csv");
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