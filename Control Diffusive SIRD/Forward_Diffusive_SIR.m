function [Rho_t, outTimes, dims, N_T, tMax, box, aLine] = Forward_Diffusive_SIR(beta,gamma,SaveFig)
% function [I_t, outTimes,I_int] = Forward_DSIR(v)
% Solves the DSIR problem on a box with periodic boundary conditions
%
%   \rho = \big( \begin{smallmatrix} S \\ I \\ R \end{smallmatrix} \big)
%
%   \frac{\partial S}{\partial t} = D \Delta S - cSI
%
%   \frac{\partial I}{\partial t} = D \Delta I + cSI - wI
%
%   \frac{\partial R}{\partial t} = D \Delta R + wI
%
%
% where for the usual SIR we would need
%
%   beta <- c, gamma <- w, D <- 0
%
if nargin == 0
  SaveFig = true;
  c = 0.479;  % Infection rate: beta
  w = 0.125;  % Recovery rate: gamma
end
if nargin == 2
    SaveFig = true;
end

    % Parameters
    D = 0.01;   % Diffusion of compartments
    c = beta;
    w = gamma;

    % number of points in y1 and y2 directions (a few 10's is normally more
    % than enough)
    N1 = 30;    % was 30
    N2 = 30;    % was 30
    N = N1 * N2;
    dims = {N1, N2, N};
    
    maskS = 1:N;
    maskI = N+1:2*N;
    maskR = 2*N+1:3*N;

    % define the geometry of the box
    geom.y1Min = -5;
    geom.y1Max =  5;
    geom.y2Min = -5;
    geom.y2Max =  5;
    geom.N = [N1;N2];

    % set up a box object using the pre-defined class
    box = BoxFourierFourier(geom);

    % create the points, differentiation matrices, integration vector, and
    % 'indices', which give the masks for the boundary, etc
    [Pts,Diff,Int,Ind] = box.ComputeAll();
    % Define some aliases
    div    = Diff.div;
    grad   = Diff.grad;
    normal = Ind.normal;
    bound  = Ind.bound;

    % uncomment this to plot the grid points
    % box.PlotGrid
    % pause

    % compute the interpolation matrix from the box points to a uniform
    % grid; this is only for plotting
    Interp = box.ComputeInterpolationMatrix(...
                                    (0:0.02:1)',(0:0.02:1)',true,true);

    % define an initial condition and retrieve initial population
    [rho_ic,N_T] = InitialCondition(Pts.y1_kv,Pts.y2_kv,0);
    % For normalised surfaces use:
    %c = c*N_T;     % Warning! c is now in [0,N_T]
    %N_T = 1.0;

    % uncomment this to plot this initial condition as a surface
    %box.plot(rho_ic(1:N));
    %pause;

    % maximum time and times for plotting
    tMax = 50.0; % Usually T = 50
    n_t  = 70;   % Was 100
    %outTimes = [0:tMax/n_t:tMax];
    % Set up a time line where we are currently at 
    ge.yMin  = 0;
    ge.yMax  = tMax;
    ge.N     = n_t;
    aLine    = SpectralLine(ge);
    outTimes = aLine.Pts.y;

    % solve the PDE, the RHS of which is given below
    % note that setting the mass matrix on the boundary allows us to solve
    % algebraic constraints (i.e., the no-fluc BC) there.

    mM          = ones([N,3]);
    mM(bound,:) = 0;
    %opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass', diag(mM(:)));
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Stats','on');
    tic
    %[~,Rho_t] = ode15s(@rhs,outTimes,rho_ic,opts);
    [~,Rho_t] = ode45(@rhs,outTimes,rho_ic,opts);
    toc

    % Get masked result for each component
    S_t = Rho_t(:, maskS)/N_T;
    I_t = Rho_t(:, maskI)/N_T;
    R_t = Rho_t(:, maskR)/N_T;
    
    if SaveFig
        % constants for plotting
        %maxRho = max(max(Rho_t));
        %minRho = min(min(Rho_t));
        %m_ic = Int*rho_ic(N+1:2*N);

        % Plot solution
        h  = figure('Position',[100,100,1200,400]);
        axis tight manual % this ensures that getframe() returns a consistent size

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
        filename = strcat( 'SIR-Forward-0_(',num2str(c),',',num2str(w),').gif' );

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
            title('$S$','Interpreter','latex')
            hold off
            % Plot I
            axes(hI)
            box.plot(It,opts);
            zlim(hI,[MinI MaxI]);   caxis(hI,[0 MaxI]);
            xlabel(''); ylabel(''); zlabel('')
            title('$I$','Interpreter','latex')
            hold off
            % Plot R
            axes(hR)
            box.plot(Rt,opts);
            zlim(hR,[MinR MaxR]);   caxis(hR,[0 MaxR]);
            xlabel(''); ylabel(''); zlabel('')
            title('$R$','Interpreter','latex')


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

    function dydt = rhs(t,rho)

        S = rho(1:N);
        I = rho(N+1:2*N);
        R = rho(2*N+1:end);

        flux = getFlux(S,I,R);

        %
        dydt = div * flux + Compartments(S,I,R);

        % no-flux BCs: gives that j.n = 0 on boundary
        dydt(bound,:) = normal * flux;

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

 function [y,N_T] = InitialCondition(y1,y2,t)
        % S initial condition is a gaussian with amplitude 7.964, variance 100/50 = 2, centred at 0
        % I initial is a mollified donut with volume (f * Int * S)
        % R initial is 0

        S = 7.964 * exp( -0.75 * ((y1-2).^2 + (y2+2).^2) );
        S = 0.7 * S + 7.964 * exp( -0.75 * ((y1+0.75).^2 + (y2-1.5).^2) ) * 0.3;
        S_T = Int * S; % Total number of people from susceptible
        
        I = exp( -1 *( (y1+0.75).^2 + (y2-1.5).^2 ));
        II = exp( -3 *( sqrt( (y1+2.5).^2 + (y2+3).^2) - 1).^2);
        I = 0.75*I + 0.25* II;
        I = S_T * I/(Int*I);
        
        % Scale
        f = 0.5; % 1e-3
        I = f * I;
        S = (1-f) * S;
        
        R = zeros([N,1]);
        
        N_T = sum([Int*S,Int*I,Int*R]);    % Compute total number of people in the population
        % Normalise
        %S = S/N_T;
        %I = I/N_T;
        y = [S;I;R];
 end

    function [y,N_T] = InitialCondition_Gaussian(y1,y2,t)
        % S initial condition is a gaussian with amplitude 7.964, variance 100/50 = 2, centred at 0
        % I initial is a mollified donut with volume (f * Int * S)
        % R initial is 0

        S = 7.964 * exp( -0.75 * (y1.^2 + y2.^2) );
        S_T = Int * S; % Total number of people from susceptible
        
        I = exp( -1 *( (y2-1.5).^2 + (y1+0.75).^2));
        I = S_T * I/(Int*I);
        
        % Normalise
        %S = S/dot(Int,S');
        f = 1e-3; % 1e-3
        I = f * I;
        S = (1-f) * S;
        
        R = zeros([N,1]);
        
        N_T = sum([Int*S,Int*I,Int*R]);    % Compute total number of people in the population
        y = [S;I;R];
    end

    function [y,N_T] = InitialCondition_Donut(y1,y2,t)
        % S initial condition is a gaussian with amplitude 7.964, variance 100/50 = 2, centred at 0
        % I initial is a mollified donut with volume (f * Int * S)
        % R initial is 0

        S = 7.964 * exp( -0.75 * (y1.^2 + y2.^2) );
        S_T = Int * S; % Total number of people from susceptible
        
        I = exp( -1 *( sqrt(y2.^2 + y1.^2) - 5).^2); % 5
        I = S_T * I/(Int*I);
        
        % Normalise
        %S = S/dot(Int,S');
        f = 1e-3; % 1e-3
        I = f * I;          % for the other I: exp( -1 *( (y2-1.5).^2 + (y1+0.75).^2));
        S = (1-f) * S;
        
        R = zeros([N,1]);
        
        N_T = sum([Int*S,Int*I,Int*R]);    % Compute total number of people in the population
        y = [S;I;R];
    end

    function [y,N_T] = InitialCondition_Scale(y1,y2,t)
        % S initial condition is a gaussian with amplitude 7.964, variance 100/50 = 2, centred at 0
        % I initial is 0.001 S
        % R initial is 0

        S = 7.964 * exp( -0.25 * (y1.^2 + y2.^2) );
        %I = 0.001 * S;
        % Normalise
        %S = S/dot(Int,S');
        f = 0.5; % 1e-3
        I = f * S;          % for the other I: exp( -1 *( (y2-1.5).^2 + (y1+0.75).^2));
        S = (1-f) * S;
        
        R = zeros([N,1]);
        
        N_T = sum([Int*S,Int*I,Int*R]);    % Compute total number of people in the population
        y = [S;I;R];
    end


end
