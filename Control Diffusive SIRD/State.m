function Rho_t = State(c,w, rho_ic, dims, box, aLine,tols)
%function [I_t, outTimes,I_int] = Forward_SIR(v)
% Solves the SIR-DDFT problem on a box with periodic boundary conditions
%
%   \rho = \big( \begin{smallmatrix} S \\ I \\ R \end{smallmatrix} \big)
%
%   \frac{\partial S}{\partial t} = D \Delta S - cSI
%
%   \frac{\partial I}{\partial t} = D \Delta I + cSI - wI
%
%   \frac{\partial R}{\partial t} = D \Delta R + wI
%
% where for the usual SIR we would need
%
%   beta <- c, gamma <- w, D <- 0
%
if nargin <= 6
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
    L      = Diff.Lap;      % for second derivative
    normal = Ind.normal;
    bound  = Ind.bound;
    % Zero matrix for blocks in Jacobian
    O  = sparse(N,N);
    oN = ones(N,1);

    % times for simulation
    outTimes = aLine.Pts.y;

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
        %[~,Rho_t] = ode15s(@rhs,outTimes,rho_ic,opts);
        [~,Rho_t] = ode45(@rhs,outTimes,rho_ic,opts);
        toc
    else
        [~,Rho_t] = ode45(@rhs,outTimes,rho_ic,opts);
    end
    
    %----------------------------------------------------------------------

    % the RHS of the PDE

    function dydt = rhs(t,rho)
        % Retrieve compartments
        S = rho(maskS);
        I = rho(maskI);
        R = rho(maskR);

        flux = getFlux(S,I,R);

        %
        dydt = div * flux + Compartments(S,I,R);

        % no-flux BCs: gives that j.n = 0 on boundary
        %dydt(bound,:) = normal * flux;

        dydt = dydt(:);         % reshape(dydt, [2*n,1]);

    end

    function f = getFlux(S,I,R)
        % the flux is just the diffusion part
        f = D * grad * [S, I, R];
    end

    function h = Compartments(S,I,R)
        cS = c * S .* I;
        cR = w * I;
        cI = cS - cR;
        h  = [-cS,cI,cR];
    end

    function J = JacobiFun(t,rho)
        % Retrieve compartments
        S = rho(maskS);
        I = rho(maskI);
        R = rho(maskR);
        
        J = [D * L - scalarOperator(c * I), scalarOperator(-c * S), O; 
             scalarOperator(c * I), D * L + scalarOperator(c * S - w), O; 
             O, scalarOperator(w * oN), D * L];
    end

    function S = scalarOperator(s)
        S = spdiags( s,0,N,N );
    end
    %scalarOperator = @(s) spdiags( s,0,N,N )
end
