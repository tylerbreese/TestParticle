classdef ShockFields

    properties
        x0 double {mustBeNumeric} = [0, 0, 0] % Shock location
        B0 double {mustBeNumeric} = [1/sqrt(2), 0, 1/sqrt(2)] % Upstream magnetic field vector
        U0 double {mustBeNumeric} = [1/sqrt(2), 0, 1/sqrt(2)] % Upstream velocity vector
        N (1,1) double {mustBeInteger, mustBePositive} = 200 % number of modes
        w_g double {mustBeNumeric, mustBePositive} = 1 % gyrofrequency
        s double {mustBeNumeric, mustBePositive} = 0.5 %variance s^2 / B^2 = number < 1
        An double {mustBeNumeric} % turbulence amplitudes
        kk double {mustBeNumeric} % k-vectors
        B double {mustBeNumeric}  % magnetic field vector at current X
        U double {mustBeNumeric}  % velocity field vector at current X
    end
    methods
        function obj = ShockFields(B0, U0, x0, N)
            if nargin > 0
                obj.B0 = B0;
                obj.U0 = U0;
                obj.x0 = x0;
                obj.N = N;
            end
        end
        
        %function [B,U] = shock_field(obj,B0,U0,x0,X,An,kk)
        function [B,U] = shock_field(obj, B0, U0, x0, X, An, kk, P, th, del)

            % for the nth time cycle update B based on location
            %if upstream of shock B0
            %if downstream of shock r*B0
            %shock taken in xz plane with downstream being the positive direction
            [r,a,b] = P.init_shock(U0, th, del);
            %obj.N = size(kk,2); 
            [dB,~] = obj.sum_turby(X,obj.N,kk,An);

            if X(1) < x0(1)
                U(1) = U0(1);
                U(2) = U0(2);
                U(3) = U0(3);
                B(1) = B0(1) + dB(1);
                B(2) = B0(2) + dB(2);
                B(3) = B0(3) + dB(3);
            else
                U(1) = U0(1)./r;
                U(2) = U0(2) + b .* U0(1);
                U(3) = U0(3) + b .* U0(1);
                B(1) = B0(1) + dB(1);
                B(2) = B0(2) + dB(2);
                B(3) = a .* ( B0(3) + dB(3) );
            end
        end

        function [An,kk] = init_turby(obj,U,w_g,s,N)
            %initialize spectrum
            %needs position, local streaming velocity and frequency for each trajectory
            %user pick Number of modes
            r_gn = vecnorm(U,2,1) ./ w_g(1);
            L = 1.496e+13; %1 AU in cm
            l1 = 5 * L; l2 = 0.5*r_gn;
            [kk,dkk] = obj.get_k_vec(l2,l1,N);

            gamma = 11/3; %3D turb`
            P0 = kk.^2 .* dkk; %make the psd array
            GkN = P0 ./ (1 + (kk.*L).^gamma);
            Gk = sum(GkN,2);

            An = 2.0 * sqrt( s^2 ./ Gk) .* sqrt(GkN);
        end

        function [K,dK] = get_k_vec(obj,l_min,l_max,N)
            %returns n by N matrix of N k vectors for n samples
            %logspace between kmin and kmax
            % excessive but i like this in a subroutine
            %N+1 so you get same modes as diff
            k_min = (2*pi)/l_max;
            k_max = (2*pi)/l_min;
            n = length(l_max);
            K = zeros(n,N+1);
            for ii = 1:n
                a = log10(k_min(ii)); b = log10(k_max(ii));
                K(ii,:) = logspace(a,b,N+1);
            end
            %get delta K
            dK = diff(K,1,2); K = K(:,2:end);
        end

        function [dB,dBn] = sum_turby(obj,X,N,kk,An)

            persistent turby

            if isempty(turby)
                [n,~,~] = size(X); %grab sample size
                beta = (2*pi).*rand(N,1); %random phase
                u = 2.*rand(n,1) - 1; % rand on -1 to 1
                %th = acos(u);
                phi = (2*pi).*rand(n,1);
                alpha = (2*pi).*rand(n,1);

                % Calculate rotation matrix components
                sin_u = sqrt(1 - u.^2);  % sqrt(1 - u^2)

                %make wavelets
                Anx =  An .* cos(alpha) .* cos(phi) .* u;
                Any =  An .* cos(alpha) .* sin(phi) .* u;
                Anz = -An .* cos(alpha) .* sin_u;
                Bnx =  An .* sin(alpha) .* sin(phi);
                Bny = -An .* sin(alpha) .* cos(phi);

                kx = kk .* cos(phi) .* sin_u;
                ky = kk .* sin(phi) .* sin_u;
                kz = kk .* u;

                % Store everything in persistent struct
                turby.kx = kx; turby.ky = ky; turby.kz = kz;
                turby.beta = beta;
                turby.Anx = Anx; turby.Bnx = Bnx;
                turby.Any = Any; turby.Bny = Bny;
                turby.Anz = Anz;
            end

            arg = (turby.kx .* X(1) + turby.ky .* X(2) + turby.kz .* X(3) + turby.beta');
            % sum
            dBn(:,:,1) = turby.Anx .* cos(arg) + turby.Bnx .* sin(arg);
            dBn(:,:,2) = turby.Any .* cos(arg) + turby.Bny .* sin(arg);
            dBn(:,:,3) = turby.Anz .* cos(arg);
            dB = squeeze(sum(dBn,2));
        end
        

    end

end