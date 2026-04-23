
folder = 'proton';
files = dir(fullfile(folder,'*h5'));
files = {files.name};
for ii = 1:length(files)
    h5files = fullfile(folder,files{ii});
    %h5disp(h5files);

    particleData = h5read(h5files,'/particle');
    m = particleData(1);
    q = particleData(2);
    Z = particleData(3);
    A = particleData(4);
    X = h5read(h5files,'/position');
    V = h5read(h5files,'/velocity');
    mask = X(:,1,2) > 0; % what is upstream of the shock?
    Vmag = vecnorm(V, 2, 2);
    Vup = squeeze(Vmag(mask,:,:));
    Vdo = squeeze(Vmag(~mask,:,:));
    
    Split_X = h5read(h5files,'/split_x');
    split_mask = Split_X(:,1) > 0;
    Split_V = h5read(h5files,'/split_p');
    Split_Vmag = vecnorm(Split_V,2,2);
    % Calculate energy for upstream and downstream particles
    En_up{ii} = 0.5 * m * Vup.^2 ./ 1.6e-6;
    En_do{ii} = 0.5 * m * Vdo.^2 ./ 1.6e-6;
    En_final{ii} = 0.5 * m * ([squeeze(Vmag(:,1,2)); Split_Vmag]).^2 ./ 1.6e-6;
    En_split{ii} = 0.5 * m * (Split_Vmag(split_mask)).^2 ./ 1.6e-6;
end

U1 = 400e5;
B1 = 0.05e-5;
bins = logspace(-4,4,100);
dE = diff(bins)./bins(2:end);
dE = dE(1); %this dE/E actually

En_up_tot = [];
En_do_tot = [];
En_total = [];
En_split_tot = [];
for jj = 1:length(files)
    En_up_tot     = cat(1,En_up_tot,En_up{jj});
    En_do_tot     = cat(1,En_do_tot,En_do{jj});
    En_total      = cat(1,En_total,En_final{jj});
    En_split_tot  = cat(1,En_split_tot,En_split{jj});
end
%%

% cutoff = [2,5,logspace(1,2,8)];
% En_cut = ((0.5 .* m .* ( U1).^2) ./ 1.6e-6) .* cutoff;
% 
% % Preallocate arrays
% numCuts = length(cutoff);
% weightedEn_split = zeros(size(En_split_tot));
% 
% for ii = 1:numCuts
%     % Logical mask for which particles exceed this energy threshold
%     mask = En_split_tot > En_cut(ii);
% 
%     % Weight for this level of splitting
%     w = 0.5^ii;
% 
%     % Apply weight only to particles that meet this condition
%     weightedEn_split(mask) = weightedEn_split(mask) + w * En_split_tot(mask);
% end
cutoff = [2, 5, logspace(1, 2, 8)];
En_cut = ((0.5 .* m .* (U1).^2) ./ 1.6e-6) .* cutoff;

numCuts = length(En_cut);

% Preallocate
weightedEn_split = zeros(size(En_split_tot));

for ii = 1:numCuts-1   % go up to numCuts-1 because we use ii and ii+1
    E_lower = En_cut(ii);
    E_upper = En_cut(ii+1);

    % Mask for energies between bounds
    mask = (En_split_tot > E_lower) & (En_split_tot <= E_upper);

    % Weight for this range (you can adjust this rule as needed)
    w = 0.5^ii;

    % Apply weighting
    weightedEn_split(mask) = weightedEn_split(mask) + w .* En_split_tot(mask);
end

% Optional: handle values above the highest cutoff
mask_high = En_split_tot > En_cut(end);
weightedEn_split(mask_high) = weightedEn_split(mask_high) + 0.5^numCuts .* En_split_tot(mask_high);

%%

plaw = bins(43:86).^(-1); plaw = plaw ./ trapz(bins(43:86),plaw);
n0 = histcounts(En_up_tot(:,3),bins);
nf = histcounts(0.5.*En_up_tot,bins);
ns = histcounts(weightedEn_split,bins);
ns = nf + ns;
n0(end+1) = 0.0; nf(end+1) = 0.0; ns(end+1) = 0.0;
figure()
hold on
grid on
xline(((0.5.*m.*(2*U1).^2)./1.6e-6),'LineWidth',2)
scatter(bins,n0./sum(n0),'o','filled','b')
%scatter(bins,nf./sum(nf),'o','filled','r')
scatter(bins,ns./sum(ns),'o','filled','g')
plot(bins(43:86),plaw,'--k','LineWidth',2)
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Kinetic Energy [MeV]')
ylabel('f(E) [Normalized Units]')
saveas(gcf,'energy.png')


%%
LECP = [4.558E+02  4.481E+00  2.147E+02  1.958E+00  1.030E+02  1.055E+00  4.924E+01  6.262E-01  1.771E+01  1.835E-01  1.261E+01  1.317E-01  5.923E+00  5.645E-02  1.358E+00  2.486E-02];
V2bins = [3.25*0.01,6*0.01, 0.1, 1.75*0.1, 3.25*0.1, 7*0.1, 1.5, 3];
jE_0 = (751.*n0) ./( (4*pi) .* (dE).*(sqrt(bins)./sqrt(2/(A*938))) );
jE_f = (751.*nf) ./( (4*pi) .* (dE).*(sqrt(bins)./sqrt(2/(A*938))) );
jE_s = (751.*ns) ./( (4*pi) .* (dE).*(sqrt(bins)./sqrt(2/(A*938))) );
pp = hex2rgb('#6F00BA');
figure()
hold on
grid on
% scatter(bins./A,jE_0.*A,'o','filled','b')
% scatter(bins./A,jE_f.*A,'o','filled','r')
% scatter(bins./A,jE_s.*A,'o','filled','g')
scatter(bins,jE_0,'o','filled','b')
%scatter(bins,jE_f,'o','filled','r')
scatter(bins,jE_s,'o','filled','g')
if contains(folder,'proton')
    scatter(V2bins,LECP(1:2:end),[],pp,'s','filled')
end
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Kinetic Energy [MeV]')
ylabel('Differential Flux J(E) [#/cm^2/s/sr/MeV/nuc]')
saveas(gcf,'flux.png')