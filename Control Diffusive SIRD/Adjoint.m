function Q_t = Adjoint(c,w, State_t, Target_t, dims, box, aLine,tols)
%
% Solves the adjoint of the SIR-DDFT problem on a box with periodic boundary conditions
%
%   q = \big( \begin{smallmatrix} q_S \\ q_I \\ q_R \end{smallmatrix} \big)
%
%   \frac{\partial q_S}{\partial t} = -D \Delta q_S + \widehat{S} - S
%                                       + cI(q_S - q_I)
%
%   \frac{\partial q_I}{\partial t} = -D \Delta q_I + \widehat{I} - I
%                                       - cI(q_S - q_I) + w(q_I - q_R)
%
%   \frac{\partial q_R}{\partial t} = -D \Delta q_R + \widehat{R} - R
%
% where for the usual SIR we would need
%
%   beta <- c, gamma <- w, D <- 0
%
if nargin <= 7
  tols = 1e-9;
end

    % Parameters
    D = 0.01;   % Diffusion of compartments

    % number of points in y1 and y2 directions (a few 10's is normally more
    % than enough)
    N = dims{3};
    
    maskS = 1:N;
    maskI = N+1:2*N;
    maskR = 2*N+1:3*N;

    % create the points, differentiation matrices, integration vector, and
    % 'indices', which give the masks for the boundary, etc
    Diff = box.ComputeDifferentiationMatrix;
    Ind  = box.ComputeIndices;
    % Define some aliases
    div    = Diff.div;
    grad   = Diff.grad;
    normal = Ind.normal;
    bound  = Ind.bound;
    L      = Diff.Lap;      % for second derivative
    % Zero matrix for blocks in Jacobian
    O  = sparse(N,N);
    oN = ones(N,1);

    % define an initial condition and retrieve initial population
    q_tc = zeros([3*N,1]);

    % times for simulation
    outTimes = aLine.Pts.y;
    revTimes = flip(outTimes);

    % solve the PDE, the RHS of which is given below
    % note that setting the mass matrix on the boundary allows us to solve
    % algebraic constraints (i.e., the no-fluc BC) there.

    %mM          = ones([N,3]);
    %mM(bound,:) = 0;
    jac   = true;
    stats = 'off';
    %opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass', diag(mM(:)));
    if jac
        opts = odeset('RelTol',tols,'AbsTol',tols,'Stats',stats,...
                        'Jacobian',@(t,rho)JacobiFun(t,rho));
    else
        opts = odeset('RelTol',tols,'AbsTol',tols,'Stats',stats);
    end
    
    timing = false;
    if timing
        tic
        [~,Q_t] = ode45(@rhs,revTimes,q_tc,opts);
        toc
    else
        [~,Q_t] = ode45(@rhs,revTimes,q_tc,opts);
    end

    % Reverse solution in time to have it back in [0,T]
    Q_t = flip(Q_t);
    
    %----------------------------------------------------------------------

    % the RHS of the PDE

    function dydt = rhs(t,q)
        % Compute interpolator at time t
        IP = aLine.ComputeInterpolationMatrixPhys(t).InterPol;
        % Interpolate target
        Target = (IP * Target_t)';
        % Interpolate state
        State = (IP * State_t)';
        
        % Retrieve each component of the target
        Ta_S = Target(maskS);   Ta_I = Target(maskI);   Ta_R = Target(maskR);
        % Retrieve each component of the state
        S = State(maskS);          I = State(maskI);       R = State(maskR);
        % Retrieve each component of the adjoint
        q_S = q(maskS);          q_I = q(maskI);         q_R = q(maskR);

        % Get flux
        flux = getFlux(q_S,q_I,q_R);
        % RHS is the diffusion term + composite terms
        dydt = div * flux + Compartments(q_S,q_I,q_R, S,I,R, Ta_S,Ta_I,Ta_R);

        % no-flux BCs: gives that j.n = 0 on boundary (redundant for periodic)
        dydt(bound,:) = normal * flux;
        dydt = dydt(:);         % reshape(dydt, [3*N,1]);

    end

    function f = getFlux(q_S,q_I,q_R)
        % the flux is just the diffusion part
        f = -D * grad * [q_S, q_I, q_R];
    end

    function h = Compartments(q_S,q_I,q_R, S,I,R, Ta_S,Ta_I,Ta_R)
        cS = Ta_S - S + c * I .* (q_S - q_I);
        cI = Ta_I - I + c * S .* (q_S - q_I) + w * (q_I - q_R);
        cR = Ta_R - R;
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
        q_S = q(maskS);          q_I = q(maskI);         q_R = q(maskR);
        
        % Assemble Jacobian
        J = [-D * L + scalarOperator(c * I), scalarOperator(-c * I), O; 
             scalarOperator(c * S), -D * L + scalarOperator(-c * S + w), scalarOperator(-w * oN); 
             O, O, -D * L];
    end

    function S = scalarOperator(s)
        S = spdiags( s,0,N,N );
    end
    %scalarOperator = @(s) spdiags( s,0,N,N )
end
