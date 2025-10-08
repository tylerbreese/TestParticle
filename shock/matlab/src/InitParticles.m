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
            Va = 30e5; % cm/s
            Vs = 70e5; % cm/s
            Ms = vecnorm(U0,2,1) ./ Vs;
            Ma = vecnorm(U0,2,1) ./ Va;
            r = ((5/3) + 1) / ((5/3) - 1 + (2/Ms^2));
            a = r * (Ma^2*cos(del)^2 - cos(th)^2) / (Ma^2*cos(del)^2 - r*cos(th)^2);
            b = ((r-1)*cos(th)*sin(th)) / (Ma^2*cos(del)^2 - r*cos(th)^2);
        end

        function V0 = sampling(obj, sample_size, U0)
            v = linspace(1000, vecnorm(U0,2,1), 100);
            u = obj.pdfrnd(linspace(-1,1,2001), ones(1,2001), sample_size);
            phi = obj.pdfrnd(linspace(0,2*pi,2001), ones(1,2001), sample_size);

            power1 = 1;
            F = v.^power1; 
            F = F ./ trapz(v,F);
            V_plasma = obj.pdfrnd(v,F,sample_size);

            vx = V_plasma .* sqrt(1 - u.^2) .* cos(phi);
            vy = V_plasma .* sqrt(1 - u.^2) .* sin(phi);
            vz = V_plasma .* u;

            V0 = [vx(:) + U0(1), vy(:) + U0(2), vz(:) + U0(3)];
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
