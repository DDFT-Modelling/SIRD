function Rho_t = State(u,v, rho_ic, dims, aLine, Conv, Diff, tols)
%
% This code does not apply the divergence to the flux explicitly, but does 
% it implicitly.
%
% Solves the SIR-DDFT problem on a box with periodic boundary conditions
%
%   \rho = \big( \begin{smallmatrix} S \\ I \\ R \end{smallmatrix} \big)
%
%   \frac{\partial S}{\partial t} =
%            \nabla \cdot \Big( D_S \nabla S - \Gamma_S S \nabla \big( K_{\text{sd}} \star (S+R) +  K_{\text{si}}
%                                \star I\big)  \Big) - cSI
%
%   \frac{\partial I}{\partial t} =
%           \nabla \cdot \Big( D_I \nabla I - \Gamma_I I \nabla \big( K_{\text{si}} \star (S+I+R) \big)  \Big)
%                + cSI - wI - mI
%
%   \frac{\partial R}{\partial t} =
%           \nabla \cdot \Big( D_R \nabla R - \Gamma_R R \nabla \big( K_{\text{sd}} \star (S+R) +  K_{\text{si}}
%                                 \star I\big)  \Big) + wI
%
%   \vec{n} \cdot \Big( D_S \nabla S - \Gamma_S S \nabla \big( u K_{sd} \star (S+R) + v K_{si} \star I\big)  \Big) = 0
%   \vec{n} \cdot \Big( D_I \nabla I - \Gamma_I I \nabla \big( v K_{si} \star (S+I+R) \big) = 0
%   \vec{n} \cdot \Big( D_R \nabla R - \Gamma_R R \nabla \big( u K_{sd} \star (S+R) + v K_{si} \star I\big)  \Big) = 0
%
% where
%
%   K_i = \exp( -s_i( x_1^2 + x_2^2) )        for all i in {sd,si}
%
% For the usual SIR we would need:
%
%   beta <- c, gamma <- w, D <- 0
%
if nargin <= 7
  tols = 1e-9;
end
if isscalar(u)
    % Scalar provided is converted into vector
    u = u * ones([aLine.N,1]);
end
if isscalar(v)
    % Scalar provided is converted into vector
    v = v * ones([aLine.N,1]);
end

    % Parameters
    % [u] Interaction strength for social distancing [p9]
    % [v] Interaction strength for self-isolation [p9] Negative->repelling
    %data.s_sd = 3;  %100.0; % Range of interaction for social distancing [p12]
    %data.s_si = 3;  %100.0; % Range of interaction for self-isolation [p12]
    % We omit the "data" structure for recursive calls
    D_S = 0.01;   % Diffusion of Susceptibles [p13]
    D_I = 0.01;   % Diffusion of Infected [p13]
    D_R = 0.01;   % Diffusion of Recovered [p13]
    c   = 0.479;  % Infection rate [p12]
    w   = 0.125;  % Recovery rate [p13]
    m   = 0.0007; % Death rate [p12]
    G_S = 1.0;    % Mobility of Susceptibles [p8]
    G_I = 1.0;    % Mobility of Infected [p8]
    G_R = 1.0;    % Mobility of Recovered [p8]

    % number of points in y1 and y2 directions (a few 10's is normally more
    % than enough)
    N = dims{3};
    
    maskS = 1:N;
    maskI = N+1:2*N;
    maskR = 2*N+1:3*N;

    % create the points, differentiation matrices, integration vector, and
    % 'indices', which give the masks for the boundary, etc
    %Diff = box.ComputeDifferentiationMatrix;
    %Ind  = box.ComputeIndices;
    % Define some aliases
    div    = Diff.div;
    grad   = Diff.grad;
    L      = Diff.Lap;      % for second derivative
    Dy_1   = Diff.Dy1;
    Dy_2   = Diff.Dy2;
    %normal = Ind.normal;
    %bound  = Ind.bound;
    % Zero matrix for blocks in Jacobian
    O  = sparse(N,N);
    oN = ones(N,1);

    % times for simulation
    outTimes = aLine.Pts.y;
    
    % compute the matrix Conv such that Conv*rho gives the
    % convolution of the kernel with rho
    %Conv = box.ComputeConvolutionMatrix(@Kernels,true);
    
    % precompute \nabla K_u (*) and \nabla K_v (*)
    Conv_u = Conv(:,:,1);
    Conv_v = Conv(:,:,2);

    % solve the PDE, the RHS of which is given below
    % note that setting the mass matrix on the boundary allows us to solve
    % algebraic constraints (i.e., the no-fluc BC) there.

    %mM          = ones([N,3]);
    %mM(bound,:) = 0;
    jac   = false;
    stats = 'off';
    %opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass', diag(mM(:)));
    if jac
        LConv_u = L * Conv(:,:,1);
        LConv_v = L * Conv(:,:,2);
        gConv_u = grad * Conv(:,:,1); % 0.046002
        gConv_v = grad * Conv(:,:,2); % 0.046002
        
        opts = odeset('RelTol',tols,'AbsTol',tols,'Stats',stats,...
                    'Jacobian',@(t,rho)JacobiFun(t,rho));
    else
        opts = odeset('RelTol',tols,'AbsTol',tols,'Stats',stats, ...
            'NonNegative',1:3*N ); 
        % This seems to be needed for D-SIR, but it is too expensive here
    end
    %disp(opts)
    
    timing = false;
    if timing
        tic
        %[~,Rho_t] = ode15s(@rhs,outTimes,rho_ic,opts);
        [~,Rho_t] = ode113(@rhs,outTimes,rho_ic,opts);
        toc
    else
        [~,Rho_t] = ode113(@rhs,outTimes,rho_ic,opts);
    end
    
    %----------------------------------------------------------------------

    % the RHS of the PDE

    function dydt = rhs(t,rho)
        %if min(rho,[],'all') < 0
        %   fprintf('Negative! %3.3f, %2.10e, %', 100*numel(rho(rho < 0.0))/numel(rho), t)
        %   rho = max(rho,0);
        %end
        %rho = max(rho,0);

        % Retrieve compartments
        S = rho(maskS);
        I = rho(maskI);
        R = rho(maskR);

        flux = getFlux(t,S,I,R);
        dydt = flux + Compartments(S,I,R);

        % no-flux BCs: gives that j.n = 0 on boundary
        %dydt(bound,:) = normal * flux;
        dydt = dydt(:);         % reshape(dydt, [2*n,1]);

    end

    function f = getFlux(t,S,I,R)
        % the flux is the sum of a diffusion part and the two-body interactions
        
        % Diffusion
        d = L * [D_S * S, D_I * I, D_R * R];
        
        % Two-body
        g = get2Body(t,S,I,R);
        f = d + g;
    end

    function g = get2Body(t,S,I,R)
        % the convolution has a vector as the kernel so we deal with the
        % two parts separately
        % K_sd is Conv(:,:,1)    &    K_si is Conv(:,:,2)
        
        % Compute interpolator at time t
        IP = aLine.ComputeInterpolationMatrixPhys(t).InterPol;
        % Interpolate u
        sd = IP * u;
        % Interpolate v
        si = IP * v;
        
        % Compute second order information
        H_u = Conv_u * (sd * (S+R) ) + Conv_v * (si * I);
        H_v = Conv_v * (si * (S + I + R) );
        
        % Mean: 1.2532
        %g1 = (Dy_1 * S).*(Dy_1 * H_u) + (Dy_2 * S).*(Dy_2 * H_u)  + (S .* (L * H_u));
        %g2 = (Dy_1 * I).*(Dy_1 * H_v) + (Dy_2 * I).*(Dy_2 * H_v)  + (I .* (L * H_v));
        %g3 = (Dy_1 * R).*(Dy_1 * H_u) + (Dy_2 * R).*(Dy_2 * H_u)  + (R .* (L * H_u));
        
        % Mean: 1.1592
        DHu1 = Dy_1 * H_u;      DHu2 = Dy_2 * H_u;  LHu = L * H_u;
        DHv1 = Dy_1 * H_v;      DHv2 = Dy_2 * H_v;  LHv = L * H_v;
        g1 = (Dy_1 * S).*DHu1 + (Dy_2 * S).*DHu2  + (S .* LHu);
        g2 = (Dy_1 * I).*DHv1 + (Dy_2 * I).*DHv2  + (I .* LHv);
        g3 = (Dy_1 * R).*DHu1 + (Dy_2 * R).*DHu2  + (R .* LHu);
        
        %g1 = [S;S] .* grad * H_u;
        %g2 = [I;I] .* gConv_v * (G_I * si * (S+I+R) );
        %g3 = [R;R] .* gConv_u * (G_R * sd * (S+R) ) + gConv_v * (G_R * si * I);
        g  = -[g1, g2, g3] .* [G_S, G_I, G_R];
    end

    function h = Compartments(S,I,R)
        cS = c * S .* I;
        cR = w * I;
        cI = cS - cR - m * I;
        h  = [-cS,cI,cR];
    end

    function J = JacobiFun(t,rho)
        % Retrieve compartments
        S = rho(maskS);
        I = rho(maskI);
        R = rho(maskR);
        
        %---- Non interaction terms ----%
        NIn_J = [D_S * L - scalarOperator(c * I), scalarOperator(-c * S), O; 
                 scalarOperator(c * I), D_I * L + scalarOperator(c * S - w - m), O; 
                 O, scalarOperator(w * oN), D_R * L];
        
        %J = NIn_J; 
        %---- Interaction terms ----%
        
        % Compute interpolator at time t
        IP = aLine.ComputeInterpolationMatrixPhys(t).InterPol;
        % Interpolate u
        sd = IP * u;
        % Interpolate v
        si = IP * v;
        
        
        % Gradient of Kernel
        %DK_u = sd * gConv_u;         % u \nabla K_u (*)     % We avoid this...
        %DK_v = si * gConv_v;         % v \nabla K_v (*)     % We avoid this...
        
        % Partial derivatives
        J_SR = div * scalarOperator2N([-sd * G_S * S; -sd * G_S * S]) * gConv_u;
        J_SS = J_SR + scalarOperator(LConv_u * (-sd * G_S * (S+R) ) + LConv_v * (-G_S * si * I ) );
        J_SI = div * scalarOperator2N([-si * G_S * S; -si * G_S * S]) * gConv_v;
        
        J_IS = div * scalarOperator2N([-si * G_I * I; -si * G_I * I]) * gConv_v;
        J_II = J_IS + scalarOperator(LConv_v * (-G_I * si * (S+I+R) ) );
        %J_IR = J_IS;
        
        J_RS = div * scalarOperator2N([-sd * G_R * R; -sd * G_R * R]) * gConv_u;
        J_RI = div * scalarOperator2N([-si * G_R * R; -si * G_R * R]) * gConv_v;
        J_RR = J_RS + scalarOperator(LConv_u * (-sd * G_R * (S+R) ) + LConv_v * (-G_R * si * I ) );
        
        % Assemble
        In_J = [ J_SS, J_SI, J_SR;
                 J_IS, J_II, J_IS;
                 J_RS, J_RI, J_RR];
        
        % Sum two parts
        J = NIn_J + In_J;
    end

    function S = scalarOperator(s)
        S = spdiags( s,0,N,N );
    end
    function S = scalarOperator2N(s)
        S = spdiags( s,0,2*N,2*N );
    end

    %scalarOperator = @(s) spdiags( s,0,N,N )
end


%     %% Profiling test: comment lines 112-119 and run this there instead
%     % Create a variable rho in the environment to run rhs:
%     % Run the function once
%     rhs(0.5,ones(3*N,1));
%     Rho_t = 0;
%     %%
