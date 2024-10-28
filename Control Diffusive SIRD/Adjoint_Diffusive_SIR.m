function Q_t = Adjoint_Diffusive_SIR(c,w,State_t, Target_t, dims, N_T, box, aLine, SaveFig)
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
if nargin <= 6
  SaveFig = true;
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
    [Pts,Diff,Int,Ind] = box.ComputeAll();
    % Define some aliases
    div    = Diff.div;
    grad   = Diff.grad;
    normal = Ind.normal;
    bound  = Ind.bound;

    % compute the interpolation matrix from the box points to a uniform
    % grid; this is only for plotting
    Interp = box.ComputeInterpolationMatrix(...
                                    (0:0.02:1)',(0:0.02:1)',true,true);

    % define an initial condition and retrieve initial population
    q_tc = TerminalCondition(Pts.y1_kv,Pts.y2_kv,0);

    % maximum time and times for plotting
    tMax     = aLine.yMax;
    outTimes = aLine.Pts.y;
    revTimes = flip(outTimes);

    % solve the PDE, the RHS of which is given below
    % note that setting the mass matrix on the boundary allows us to solve
    % algebraic constraints (i.e., the no-fluc BC) there.

    mM          = ones([N,3]);
    mM(bound,:) = 0;
    %opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass', diag(mM(:)));
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Stats','on');
    tic
    %[~,Q_t] = ode15s(@rhs,revTimes,q_tc,opts);
    [~,Q_t] = ode45(@rhs,revTimes,q_tc,opts);
    toc

    % Reverse solution in time to have it back in [0,T]
    Q_t = flip(Q_t);
    % Get masked result for each component
    S_t = Q_t(:, maskS)/N_T;
    I_t = Q_t(:, maskI)/N_T;
    R_t = Q_t(:, maskR)/N_T;
    
    if SaveFig
        % Plot solution
        h  = figure('Position',[100,100,1200,400]);
        axis tight manual % this ensures that getframe() returns a consistent size

        % constants for plotting
        MaxS = max(S_t,[],'all');   MinS = min(S_t,[],'all');
        MaxI = max(I_t,[],'all');   MinI = min(I_t,[],'all');
        MaxR = max(R_t,[],'all');   MinR = min(R_t,[],'all');
        % Add a little extra of height with the same magnitude
        MaxS = MaxS + 10^floor( log10(abs(MaxS)) );
        MaxI = MaxI + 10^floor( log10(abs(MaxI)) );
        MaxR = MaxR + 10^floor( log10(abs(MaxR)) );

        hS = subplot(1,3,1);  zlim(hS,[MinS MaxS]);  caxis(hS,[0 MaxS]);
        hI = subplot(1,3,2);  zlim(hI,[MinI MaxI]);  caxis(hI,[0 MaxI]);
        hR = subplot(1,3,3);  zlim(hR,[MinR MaxR]);  caxis(hR,[0 MaxR]);

        set(gca,'Color','none');   set(gcf, 'color', 'white');

        opts = {};  % plotting options - the default are ok for now
        filename = strcat( 'SIR-Adjoint-0_(',num2str(c),',',num2str(w),').gif' );

        for iTime = 1:length(outTimes)

            % extract the solution at this time
            St = S_t(iTime,:)';
            It = I_t(iTime,:)';
            Rt = R_t(iTime,:)';
            I_int(iTime) = Int*It;%sum(It*Int,'all');
            S_int(iTime) = Int*St;
            R_int(iTime) = Int*Rt;
            %rhot = rhot(n+1:2*n);  % If line 103 is commented

            % plot with built-in class function
            hold off
            % Plot S
            axes(hS)
            box.plot(St,opts);
            zlim(hS,[MinS MaxS]);   caxis(hS,[0 MaxS]);
            xlabel(''); ylabel(''); zlabel('')
            title('$q_S$','Interpreter','latex')
            hold off
            % Plot I
            axes(hI)
            box.plot(It,opts);
            zlim(hI,[MinI MaxI]);   caxis(hI,[0 MaxI]);
            xlabel(''); ylabel(''); zlabel('')
            title('$q_I$','Interpreter','latex')
            hold off
            % Plot R
            axes(hR)
            box.plot(Rt,opts);
            zlim(hR,[MinR MaxR]);   caxis(hR,[0 MaxR]);
            xlabel(''); ylabel(''); zlabel('')
            title('$q_R$','Interpreter','latex')


            % set title with some info
            sgtitle(['$t$ = ' num2str(outTimes(iTime))],'Interpreter','latex','FontSize',20)

            % fix axes
    %         colorbar

            %pause(0.05)

            % Capture the plot as an image 
            frame = getframe(h); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 

            % Write to the GIF File 
            if iTime == 1 
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
            end 

        end
    end
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

        % no-flux BCs: gives that j.n = 0 on boundary
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

    function y = TerminalCondition(y1,y2,t)
        y = zeros([3*N,1]);
    end

end
