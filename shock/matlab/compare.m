    %%grab multiple datasets and plot them to compare
    
    Folders = [
        "1000 cycles/proton"
        "1000 cycles/helium"
        "10000 cycles/helium"
        "helium"
        ];
    for kk = 1:length(Folders)
        folder = Folders(kk);
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
        for jj = 1:10
            En_up_tot     = cat(1,En_up_tot,En_up{jj});
            En_do_tot     = cat(1,En_do_tot,En_do{jj});
            En_total      = cat(1,En_total,En_final{jj});
            En_split_tot  = cat(1,En_split_tot,En_split{jj});
        end
        % weight split particles
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
        % calc distributions
        plaw = bins(43:86).^(-1); plaw = plaw ./ trapz(bins(43:86),plaw);
        n0 = histcounts(En_up_tot(:,3),bins);
        nf = histcounts(0.5.*En_up_tot,bins);
        ns = histcounts(weightedEn_split,bins);
        ns = nf + ns;
        n0(end+1) = 0.0; nf(end+1) = 0.0; ns(end+1) = 0.0;
        jE_0 = (n0) ./( (4*pi) .* (dE).*(sqrt(bins)./sqrt(2/(A*938))) );
        %jE_f = (nf) ./( (4*pi) .* (dE).*(sqrt(bins)./sqrt(2/(A*938))) );
        jE_s = (ns) ./( (4*pi) .* (dE).*(sqrt(bins)./sqrt(2/(A*938))) );
        %if per nucleon Ebins = Ebins/nuc, J' = nuc*J
    
        out{kk} = {[bins;jE_0;jE_s]};
    
    end
    
%% external data

nH = 0.115; %proton pUI density
nHe = 0.01; %helium density ratio? what is good value find the paper
%xs = 5e-15./sqrt(1e7.*bins);
%ENA calculation J_ENA = J_PUI * (dR/lambda)
a = [4.15, 0.531, 67.3];
xs = @(E) 1e-16 .* (a(1) - a(2) .* log(E)).^2 .* (1 - exp(-a(3)./(E))).^(4.5);
lambda = 1 ./ (xs(1000.*bins).*nH);
DR = 35 .* 1.496e13;
LECP_3day = [%242;243;244;
4.558E+02  4.481E+00  2.147E+02  1.958E+00  1.030E+02  1.055E+00  4.924E+01  6.262E-01  1.771E+01  1.835E-01  1.261E+01  1.317E-01  5.923E+00  5.645E-02  1.358E+00  2.486E-02
4.529E+02  4.003E+00  2.175E+02  1.766E+00  1.139E+02  9.910E-01  5.855E+01  6.094E-01  2.183E+01  1.818E-01  1.314E+01  1.200E-01  5.788E+00  4.980E-02  1.319E+00  2.186E-02
4.713E+02  5.047E+00  2.207E+02  2.198E+00  1.112E+02  1.209E+00  5.749E+01  7.450E-01  2.106E+01  2.203E-01  1.223E+01  1.428E-01  5.667E+00  6.079E-02  1.366E+00  2.745E-02];
LECP = sum(LECP_3day,1)./3;
V2bins = [3.25*0.01,6*0.01, 0.1, 1.75*0.1, 3.25*0.1, 7*0.1, 1.5, 3];
%V2lambda = 1 ./ (nH .* 5e-15./sqrt(1e7.*V2bins));
V2lambda = 1 ./ (nH .* xs(1000.*V2bins));
V2ENA = LECP(1:2:16).* (DR ./ V2lambda);
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

%% plotting
% some colors
dg = hex2rgb('#006400');
pp = hex2rgb('#6F00BA');
or = hex2rgb('#FF9800');
bb = hex2rgb('#4472C4');

nuc = A;
nH = 1;
nHe = 0.1.*nH;
dist1 = [out{1}{1}];
dist2 = [out{2}{1}];
dist3 = [out{3}{1}];
dist4 = [out{4}{1}];

figure()
hold on
grid on

scatter(dist2(1,:),dist2(2,:).*nHe,'o','filled','b')
scatter(dist2(1,:),dist2(3,:).*nHe,'o','filled','r')
scatter(dist3(1,:),10.*dist3(3,:).*nHe,'o','filled','g')

scatter(dist4(1,:),10.*dist4(3,:).*nHe,'o','filled','k')

scatter(dist1(1,:),dist1(2,:).*nH,'x','k') 
scatter(dist1(1,:),dist1(3,:).*nH,'x','m')

% scatter(dist2(1,:)./nuc,0.1.*dist2(2,:).*nuc,'o','filled','b')
% scatter(dist2(1,:)./nuc,0.1.*dist2(3,:).*nuc,'o','filled','r')
% scatter(dist3(1,:)./nuc,0.1.*10.*dist3(3,:).*nuc,'o','filled','g')
    scatter(V2bins,LECP(1:2:end),[],pp,'s','filled')

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Kinetic Energy [MeV]')
ylabel('Differential Flux J(E) [#/cm^2/s/sr/MeV]')

legend('He+, Initial','He+, Final (1000\Omega_g^{-1})','He+, Final (10000\Omega_g^{-1})','He+, Final (100000\Omega_g^{-1})', ...
    'H+, Initial (1000\Omega_g^{-1})','H+, Final (1000\Omega_g^{-1})')
