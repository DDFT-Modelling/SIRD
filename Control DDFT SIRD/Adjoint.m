function Q_t = Adjoint(u,v, State_t, dims, aLine, Conv, Diff, tols)
%
% Solves the adjoint of the SIR-DDFT problem on a box with periodic boundary conditions
%
%   q = \big( \begin{smallmatrix} q_S \\ q_I \\ q_R \end{smallmatrix} \big)
%
%   \frac{\partial q_S}{\partial t} = -\nabla \cdot ( D_S \nabla q_S)  
%      + cI(q_S - q_I)
%      + u \nabla K_u (*) ( \Gamma_S S \nabla q_S + \Gamma_R R \nabla q_R )
%      + v \nabla K_v (*) ( \Gamma_I I \nabla q_I)
%      - \Gamma_S \nabla q_S \cdot \nabla ( u K_u \star (S+R) + v K_v \star I )
% 
%   \frac{\partial q_I}{\partial t} = -\nabla \cdot (D_I \nabla q_I) - I
%      + cS(q_S - q_I) + w(q_I - q_R) + mq_I
%      + v \nabla K_v (*) ( \Gamma_S S \nabla q_S + \Gamma_I I \nabla q_I + \Gamma_R R \nabla q_R )
%      - \Gamma_I \nabla q_I \cdot \nabla ( v K_v \star (S+I+R) )
% 
%   \frac{\partial q_R}{\partial t} = -\nabla \cdot ( D_R \nabla q_R) 
%      + u \nabla K_u (*) ( \Gamma_S S \nabla q_S + \Gamma_R R \nabla q_R )
%      + v \nabla K_v (*) ( \Gamma_I I \nabla q_I)
%      - \Gamma_R \nabla q_R \cdot \nabla ( u K_u \star (S+R) + v K_v \star I )
%
%   \vec{n} \cdot \big( \nabla q_S \big) = 0
%   \vec{n} \cdot \big( \nabla q_I \big) = 0
%   \vec{n} \cdot \big( \nabla q_R \big) = 0
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
if length(u) == 1
    % Scalar provided is converted into vector
    u = u * ones([aLine.N,1]);
end
if length(v) == 1
    % Scalar provided is converted into vector
    v = v * ones([aLine.N,1]);
end

    % Parameters
    % [u] Interaction strength for social distancing [p9]
    % [v] Interaction strength for self-isolation [p9] Negative->repelling
    %data.s_sd = 3.0;  %100.0; % Range of interaction for social distancing [p12]
    %data.s_si = 3.0;  %100.0; % Range of interaction for self-isolation [p12]
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
    %Diff = box.ComputeDifferentiationMatrix;        % [Profiler] expensive
    %Ind  = box.ComputeIndices;
    % Define some aliases
    div    = Diff.div;
    grad   = Diff.grad;
    L      = Diff.Lap;      % for second derivative
    %normal = Ind.normal;
    %bound  = Ind.bound;
    % Zero matrix for blocks in Jacobian
    O  = sparse(N,N);
    oN = ones(N,1);

    % define an initial condition and retrieve initial population
    q_tc = zeros([3*N,1]);

    % times for simulation
    outTimes = aLine.Pts.y;
    revTimes = flip(outTimes);
    
    % compute the matrix Conv such that Conv*rho gives the
    % convolution of the kernel with rho
    %Conv = box.ComputeConvolutionMatrix(@Kernels,true);
    % precompute \nabla K_u (*) and \nabla K_v (*)
    gConv_u = grad * Conv(:,:,1); %786 * 0.046002       [Profiler] expensive
    gConv_v = grad * Conv(:,:,2); %786 * 0.046002       (would it be a good idea to pass all these as arguments in a structure?)
    gCu_x = gConv_u(1:N,:);       % 0.029149
    gCu_y = gConv_u(N+1:2*N,:);   % 0.013642
    gCv_x = gConv_v(1:N,:);       % 0.029149
    gCv_y = gConv_v(N+1:2*N,:);   % 0.013642
    

    % solve the PDE, the RHS of which is given below
    % note that setting the mass matrix on the boundary allows us to solve
    % algebraic constraints (i.e., the no-fluc BC) there.

    %mM          = ones([N,3]);
    %mM(bound,:) = 0;
    jac   = false;
    stats = 'off';
    %opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass', diag(mM(:)));
    if jac
        grad_x = Diff.Dy1;      %grad(1:N,:);
        grad_y = Diff.Dy2;      %grad(N+1:end,:);
        opts = odeset('RelTol',tols,'AbsTol',tols,'Stats',stats,...
                        'Jacobian', @(t,rho)JacobiFun(t,rho));
    else
        opts = odeset('RelTol',tols,'AbsTol',tols,'Stats',stats);
    end
    
    timing = false;
    if timing
        tic
        [~,Q_t] = ode113(@rhs,revTimes,q_tc,opts);
        toc
    else
        [~,Q_t] = ode113(@rhs,revTimes,q_tc,opts);
    end

    % Reverse solution in time to have it back in [0,T]
    Q_t = flip(Q_t);
    
    %----------------------------------------------------------------------

    % the RHS of the PDE

    function dydt = rhs(t,q)
        % Tests of time for 40x40
        
        % Compute interpolator at time t
        IP = aLine.ComputeInterpolationMatrixPhys(t).InterPol;
        % Interpolate state
        State = (IP * State_t)';
        % Interpolate u and v
        sd = IP * u;
        si = IP * v;
        
        % Retrieve each component of the state
        S = State(maskS);          I = State(maskI);       R = State(maskR);
        % Retrieve each component of the adjoint
        q_S = q(maskS);          q_I = q(maskI);         q_R = q(maskR);
        % Matricise vector to obtain derivatives (recycle variables)
        Q = reshape(q, [N,3]);      % 0.006711
        gQ = grad * Q;              % 0.006182
        
        % Get flux
        flux = getFlux(gQ);
        % Advection
        adv = Advection(gQ, S,I,R, sd,si);
        
        % RHS is the diffusion term + composite terms
        dydt = div * flux + Compartments(q_S,q_I,q_R, S,I,R) + adv;

        % no-flux BCs: gives that j.n = 0 on boundary
        %dydt(bound,:) = normal * flux;
        dydt = dydt(:);         % reshape(dydt, [3*N,1]);

    end

    function f = getFlux(gQ)
        % the flux is just the diffusion part
        %f = -grad * [D_S * q_S, D_I * q_I, D_R * q_R];     % 
        f = gQ .* [-D_S -D_I -D_R];                         % 0.00563
    end

    % [Profiler] expensive
    function g = Advection(gQ, S,I,R, sd, si)%(q_S,q_I,q_R, gQ, S,I,R, sd,si)
        % Gradient of Kernel
        %DK_u = sd * gConv_u;         % u \nabla K_u (*)     % We avoid this...
        %DK_v = si * gConv_v;         % v \nabla K_v (*)     % We avoid this...
        
        
        
        % State times gradients of adjoints
        G_S_Dq = G_S * [S;S] .* gQ(:,1);     % \Gamma_S S \nabla q_S
        G_I_Dq = G_I * [I;I] .* gQ(:,2);     % \Gamma_I I \nabla q_I
        G_R_Dq = G_R * [R;R] .* gQ(:,3);     % \Gamma_R R \nabla q_R
        
        % Partial sums
        G_SR_Dq  = G_S_Dq + G_R_Dq;
        G_SIR_Dq = G_I_Dq + G_SR_Dq;
        
        % Products to be recycled
        %A_SR_u = sd * (gConv_u(1:N,:) * G_SR_Dq(1:N) + gConv_u(N+1:end,:) * G_SR_Dq(N+1:end));      % 0.016084
        A_SR_u = sd * (gCu_x * G_SR_Dq(1:N) + gCu_y * G_SR_Dq(N+1:2*N));                             % 0.005982
        A_I_v  = si * (gCv_x *  G_I_Dq(1:N) + gCv_y *  G_I_Dq(N+1:2*N));                             % 0.010024
        
        % ------------------
        % Advection in q_S: Express dot product as Hadamard -> sum
        %AS = -(grad * G_S * q_S) .* ( grad * (sd*Conv(:,:,1) * (S+R) + si*Conv(:,:,2) * I ) );      % 0.017016
        aux_AS = gConv_u * (sd * (S+R)) + gConv_v * (si * I);                                        % 0.006742
        AS = -G_S * gQ(:,1) .* aux_AS;                                                               % 0.009166 -> 0.005176
        % Sum components of Hadamard product  and add iteractions from
        % social distancing and self-isolation
        aS = AS(1:N) + AS(N+1:2*N)  +  A_SR_u + A_I_v;                                                % 0.007088
        
        % ------------------
        % Advection in q_I:
        AI = -G_I * gQ(:,2) .* ( gConv_v * (si * (S+I+R) ) );                                        % 0.009328
        % Sum components and add interactions
        aI = AI(1:N) + AI(N+1:2*N)  +  si * (gCv_x * G_SIR_Dq(1:N) + gCv_y * G_SIR_Dq(N+1:2*N));     % 0.009385
        
        % ------------------
        % Advection in q_R:
        AR = -G_R * gQ(:,3) .* aux_AS;                                                               % 0.005176
        % Sum components and add interactions
        aR = AR(1:N) + AR(N+1:2*N)  +  A_SR_u + A_I_v;
        
        % Assemble
        g = [aS,aI,aR];
    end

    function h = Compartments(q_S,q_I,q_R, S,I,R)
        cS = c * I .* (q_S - q_I);
        cI = - I + c * S .* (q_S - q_I) + w * (q_I - q_R) + m * q_I;
        cR = q_tc(maskR);
        h  = [cS,cI,cR];
    end

    function J = JacobiFun(t,q)
        % Compute interpolator at time t
        IP = aLine.ComputeInterpolationMatrixPhys(t).InterPol;
        % Interpolate state
        State = (IP * State_t)';
        
        % Retrieve each component of the state
        S = State(maskS);          I = State(maskI);       R = State(maskR);
        % Retrieve each component of the adjoint
        %q_S = q(maskS);          q_I = q(maskI);         q_R = q(maskR);
        
        %---- Non interaction terms ----%
        NIn_J = [-D_S * L + scalarOperator(c * I), scalarOperator(-c * I), O; 
                    scalarOperator(c * S), -D_I * L + scalarOperator(-c * S + (w+m)), scalarOperator(-w * oN); 
                    O, O, -D_R * L];
                
        %J = NIn_J; 
        
        %---- Interaction terms ----%

        % Interpolate u
        sd = IP * u;
        % Interpolate v
        si = IP * v;
        
        
        % Gradient of Kernel
        %DK_u = sd * gConv_u;         % u \nabla K_u (*)     % We avoid this...
        %DK_v = si * gConv_v;         % v \nabla K_v (*)     % We avoid this...
        
        % Partial derivatives
        J_SS = gCu_x * (sd * G_S * S) .* grad_x + gCu_y * (sd * G_S * S) .* grad_y + ...
                    scalarOperator2Dif(gConv_u * (-sd * G_S * (S+R) ) + gConv_v * (-G_S * si * I ) ) * grad;
        %           scalarOperator(LConv_u * (-sd * G_S * (S+R) ) + LConv_v * (-G_S * si * I ) );
        %                           <Ah, g> = <h, A'g> = <h, g'A>
        J_SI = gCv_x * (si * G_I * I) .* grad_x + gCv_y * (si * G_I * I) .* grad_y;
        J_SR = gCu_x * (sd * G_R * R) .* grad_x + gCu_y * (sd * G_R * R) .* grad_y;
                    
        J_IS = gCv_x * (si * G_S * S) .* grad_x + gCv_y * (si * G_S * S) .* grad_y;
        J_II = gCv_x * (si * G_I * I) .* grad_x + gCv_y * (si * G_I * I) .* grad_y + ...
                    scalarOperator2Dif(gConv_v * (-si * G_I * (S+I+R) ) ) * grad;
        %           scalarOperator(LConv_v * (-G_I * si * (S+I+R) ) );
        J_IR = gCv_x * (si * G_R * R) .* grad(1:N,:) + gCv_y * (si * G_R * R) .* grad_y;
        
        
        J_RS = gCu_x * (sd * G_S * S) .* grad_x + gCu_y * (sd * G_S * S) .* grad_y;
        J_RI = gCv_x * (si * G_I * I) .* grad_x + gCv_y * (si * G_I * I) .* grad_y;
        J_RR = gCu_x * (sd * G_R * R) .* grad_x + gCu_y * (sd * G_R * R) .* grad_y + ...
                    scalarOperator2Dif(gConv_u * (-sd * G_R * (S+R) ) + gConv_v * (-G_R * si * I ) ) * grad;
        %           scalarOperator(LConv_u * (-sd * G_R * (S+R) ) + LConv_v * (-G_R * si * I ) );
        
        % Assemble
        In_J = [ J_SS, J_SI, J_SR;
                 J_IS, J_II, J_IR;
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
    function S = scalarOperator2Dif(s)
        S = spdiags( [s(1:N), s(N+1:2*N)], [0,N], N, 2*N );
    end
    %scalarOperator = @(s) spdiags( s,0,N,N )
    %scalarOperator2Dif = @(s) spdiags( [s(1:N), s(N+1:2*N)], [0,N], N, 2*N );
    
end



    %% Profiling test: comment lines 123-130 and run this there instead
    % Create a variable q in the environment to run rhs:
%     IP = aLine.ComputeInterpolationMatrixPhys(0.5).InterPol;
%     q = (IP * State_t)';
%     % Run the function once
%     rhs(0.5,q);
%     Q_t = 0;
    %%