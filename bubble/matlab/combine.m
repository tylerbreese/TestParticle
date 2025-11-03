
folder = 'multi';
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

U1 = 2000e5;
B1 = 0.5;
bins = logspace(-2,1,100);
En_up_tot = [];
En_do_tot = [];
En_total = [];
for jj = 1:10
    En_up_tot     = cat(1,En_up_tot,En_up{jj});
    En_do_tot     = cat(1,En_do_tot,En_do{jj});
    En_total      = cat(1,En_total,En_final{jj});
end
%%

plaw = bins(43:86).^(-1); plaw = plaw ./ trapz(bins(43:86),plaw);
n0 = histcounts(En_up_tot(:,3),bins);
nf = histcounts(En_up_tot,bins);
% ns = histcounts(weightedEn_split,bins);
% ns = nf + ns;
n0(end+1) = 0.0; nf(end+1) = 0.0; %ns(end+1) = 0.0;
figure()
hold on
grid on
scatter(bins,n0./sum(n0),'o','filled','b')
scatter(bins,nf./sum(nf),'o','filled','r')
%scatter(bins,ns./sum(ns),'o','filled','g')
plot(bins(43:86),plaw,'--k','LineWidth',2)
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Kinetic Energy [MeV]')
ylabel('f(E) [Normalized Units]')
saveas(gcf,'energy.png')


%%
jE_0 = (n0) ./( (4*pi) .* (0.03).*(sqrt(bins)./sqrt(2/(A*938))) );
jE_f = (nf) ./( (4*pi) .* (0.03).*(sqrt(bins)./sqrt(2/(A*938))) );
%jE_s = (ns) ./( (4*pi) .* (0.03).*(sqrt(bins)./sqrt(2/(A*938))) );

figure()
hold on
grid on
scatter(bins,jE_0,'o','filled','b')
scatter(bins,jE_f,'o','filled','r')
%scatter(bins,jE_s,'o','filled','g')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Kinetic Energy [MeV]')
ylabel('Differential Flux J(E) [#/cm^2/s/sr/MeV]')
saveas(gcf,'flux.png')