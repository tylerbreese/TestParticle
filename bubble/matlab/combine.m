
folder = 'plane';
files = dir(fullfile(folder,'*h5'));
files = {files.name};
c = 3e10;
mc2 = 0.512 .* 1.6e-6;
for ii = 1:length(files)
    h5files = fullfile(folder,files{ii});
    %h5disp(h5files);

    particleData = h5read(h5files,'/particle');
    m = particleData(1);
    q = particleData(2);
    Z = particleData(3);
    A = particleData(4);
    X = h5read(h5files,'/position');
    p = h5read(h5files,'/velocity');
    mask = X(:,1,2) > 0; % what is downstream of the shock?
    pmag = vecnorm(p, 2, 2);
    pup = squeeze(pmag(mask,:,:));
    pdo = squeeze(pmag(~mask,:,:));
   
    En_up{ii}    =  sqrt( (pup.^2 .* c^2) + mc2^2 )./ 1.6e-6;
    En_do{ii}    =  sqrt( (pdo.^2 .* c^2) + mc2^2 )./ 1.6e-6;
    En_final{ii} =  sqrt( (squeeze(pmag).^2 .* c^2) + mc2^2 )./ 1.6e-6;

end
disp(files{ii})
U1 = 2000e5;
B1 = 0.5;

bins = logspace(-2,1,100);

En_up_tot = [];
En_do_tot = [];
En_total = [];
for jj = 1:length(files)
    En_up_tot     = cat(1,En_up_tot,En_up{jj});
    En_do_tot     = cat(1,En_do_tot,En_do{jj});
    En_total      = cat(1,En_total,En_final{jj});
end
%%

plaw = exp(-0.575.*bins - 0.055.*bins.^2); %plaw = plaw ./ sum(plaw);
n0 = histcounts(En_total(:,3),bins);
nf = histcounts(En_total(:,2),bins);
% ns = histcounts(weightedEn_split,bins);
% ns = nf + ns;
n0(end+1) = 0.0; nf(end+1) = 0.0; %ns(end+1) = 0.0;
figure()
hold on
grid on
scatter(bins,n0./sum(n0),'o','filled','b')
scatter(bins,nf./sum(nf),'o','filled','r')
%scatter(bins,ns./sum(ns),'o','filled','g')
%plot(bins,plaw,'--k','LineWidth',2)
ax = gca;
%ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Kinetic Energy [MeV]')
ylabel('f(E) [Normalized Units]')
saveas(gcf,'energy.png')
%%
figure()
hold on
grid on
scatter(squeeze(X(:,1,2))./1e5, En_total(1:1e4,2), [],'r','x' )
scatter(squeeze(X(:,1,3))./1e5, En_total(1:1e4,3), [],'b','filled','o' )
xlabel('Radial Position (km)')
ylabel('Kinetic Energy (MeV)')
saveas(gcf,'space.png')
%%
figure()
polarscatter(En_total(1:1e4,3),squeeze(X(:,1,3))./1e5 )
hold on
polarscatter(En_total(1:1e4,2),squeeze(X(:,1,2))./1e5 )