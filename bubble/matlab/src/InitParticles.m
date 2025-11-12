classdef InitParticles
    properties
        U0 double {mustBeNumeric} = [1/sqrt(2), 0, 1/sqrt(2)]  % Upstream velocity
        sample_size (1,1) double {mustBeInteger, mustBePositive} = 1000  % Number of particles
        del (1,1) double {mustBeNumeric} = 0.0
        th  (1,1) double {mustBeNumeric} = pi/2
    end

    methods
        function obj = InitParticles(U0, th, del, sample_size)
            if nargin > 0
                obj.U0 = U0;
                obj.th = th;
                obj.del = del;
                obj.sample_size = sample_size;
            end
        end

        function [r,a,b] = init_shock(obj,U0, th, del)
            Va = 0.6e5; % cm/s
            Vs = 300e5; % cm/s
            Ms = vecnorm(U0,2,1) ./ Vs;
            Ma = vecnorm(U0,2,1) ./ Va;
            %r = ((1.4) + 1) / ((1.4) - 1 + (2/Ms^2));
            r = 6;
            a = r * (Ma^2*cos(del)^2 - cos(th)^2) / (Ma^2*cos(del)^2 - r*cos(th)^2);
            b = ((r-1)*cos(th)*sin(th)) / (Ma^2*cos(del)^2 - r*cos(th)^2);
        end

        function p0 = sampling(obj,sample_size)
            c = 3e10;
            m = 9.11e-28;
            %beta energy dist, assume U235
            W = logspace(log10(0.1),log10(10.1),sample_size);
            F = exp(-0.575.*W - 0.055.*W.^2);
            E_plasma = datasample(W,sample_size,'Weights',F);
            E_plasma = 1.6e-6 .* E_plasma; %convert to ergs
            g = sqrt(1 + (E_plasma) ./ (m*c*c));
            p_plasma = (m*c) .* (g.^2 -1).^(0.5);
            %spherical particle distribution
            u = obj.pdfrnd(-1:0.001:1,ones(size(-1:0.001:1)),sample_size);
            phi_max = 2*pi; phi_min = 0;
            div = abs(phi_max-phi_min)/length(u); %force square array of samples
            phi = obj.pdfrnd(phi_min:div:phi_max,ones(size(phi_min:div:phi_max)),sample_size);
            p0 = [
                ( p_plasma .* ( (1 - u.^2).^0.5 .* cos(phi) )' )
                ( p_plasma .* ( (1 - u.^2).^0.5 .* sin(phi) )' )
                ( p_plasma .* u' )
                ];
        end

        function X = pdfrnd(obj,x, px, sample_size)
            cdf = cumsum(px);
            cdf = cdf / sum(cdf);
            cdf = cdf / max(cdf);
            rnd = rand(sample_size, 1);
            X = interp1(cdf, x, rnd, 'linear', 0);
        end
    end
end
