function MI = MI_DiffusionExact(opts)
   
    AddPaths();
    
    close all;  

    if(nargin==0)
        opts = struct;
    end
    
    if(isfield(opts,'BC'))
        BC = opts.BC;
    else
        BC = 'Prescribed';
    end
    
    if(isfield(opts,'endpoints'))
        endpoints = opts.endpoints;
    else
        endpoints = [-1;-0.5;0.5;1];
    end

    nIntervals = length(endpoints)-1;
    
    if(isfield(opts,'N'))
        N = opts.N;
    else
        N = 10;
    end

    if(length(N)==1)
        N = repmat(N,nIntervals,1);
    end
    
    %----------------------------------------------------------------------
    % Initialisation
    %----------------------------------------------------------------------
        
    intervals(nIntervals).interval = []; % preallocate
    
    for iInterval = 1:nIntervals
        intervals(iInterval).interval = 'SpectralLine';
        geom.N    = N(iInterval);
        geom.yMin = endpoints(iInterval);
        geom.yMax = endpoints(iInterval+1);
        intervals(iInterval).geom = geom;
    end
    
    MI = MultiInterval(intervals);
        
    BCs = struct;
    
    for iInterval = 1:nIntervals-1
        BCs(iInterval,iInterval+1).function = 'BCmatch';
    end
    
    Pts = MI.Pts;
    Diff = MI.Diff;
    Dy = Diff.Dy;
    DDy = Diff.DDy;
    Int = MI.Int;
    Ind = MI.Ind;
    
    y = Pts.y;
    yMin = Pts.yMin;
    yMax = Pts.yMax;
    bound = Ind.bound;
    intersections = Ind.intersections;
        
    %----------------------------------------------------------------------
    % Find initial condition
    %----------------------------------------------------------------------

    t0 = 0.1;

    switch BC
        case 'Prescribed'
            rho_ic = ExactSolutionPrescribed(y,t0);
        case 'Constant'
            rho_ic = ExactSolutionConstant(y,t0);
        case 'NoFlux'
            rho_ic = ExactSolutionNoFlux(y,t0);
    end
    
    %---------------------------------------------------------------------%
    % Solve the PDE                                                       %
    %---------------------------------------------------------------------%
        
    outTimes = 0.1:0.05:1;
    
    mM                = ones(size(y));
    mM(bound)         = 0;    
    mM(intersections) = 0;
    opts = odeset('RelTol',10^-12,'AbsTol',10^-12,'Mass',diag(mM));
    
    tic
    switch BC
        case 'Prescribed'
            [outTimes,rho_t] = ode15s(@rhsPrescribed,outTimes,rho_ic,opts);
        case 'Constant'
            [outTimes,rho_t] = ode15s(@rhsConstant,outTimes,rho_ic,opts);
        case 'NoFlux'
            [outTimes,rho_t] = ode15s(@rhsNoFlux,outTimes,rho_ic,opts);
    end
    odeTime = toc;
    
    disp(['Dynamics computation time = ' num2str(odeTime) 's']);
    
    nT = length(outTimes);
    
    m = zeros(nT,1);
    err = zeros(nT,1);

    rhoMax = max(max(abs(rho_t)));
    
    for iT = 1:nT
        rhot = rho_t(iT,:)';
        
        m(iT) = Int*rhot;
        
        switch BC
            case 'Prescribed'
                rhoEt = ExactSolutionPrescribed(y,outTimes(iT));
            case 'Constant'
                rhoEt = ExactSolutionConstant(y,outTimes(iT));
            case 'NoFlux'
                rhoEt = ExactSolutionNoFlux(y,outTimes(iT));
        end
        
        err(iT) = max(abs(rhoEt-rhot)./abs(rhoEt));
   
        hold off
        MI.Plot(rhot);
        
        ylim([0,rhoMax]);
        
        pause(0.1);
        
        
    end
            
    disp(['Max relative error = ' num2str(max(err))]); 
    
    switch BC
        case 'NoFlux'
            m_ic = Int*rho_ic;
            merr = max(abs(m - m_ic)/abs(m_ic));
            disp(['Max mass error = ' num2str(merr)]);
    end
    
    %----------------------------------------------------------------------
    % RHS of ODE
    %----------------------------------------------------------------------

    function drhodt = rhsPrescribed(t,rho)
                                
        drhodt  = DDy*rho;
        flux    = -Dy*rho;        

        rhoExactBound = ExactSolutionPrescribed(y(bound),t);
        drhodt(bound) = rho(bound) - rhoExactBound;
                
        data.rho  = rho;
        data.flux = flux;
        drhodt    = MI.ApplyIntersectionBCs(drhodt,data,BCs);
        
    end

    function drhodt = rhsConstant(t,rho)
                                
        drhodt  = DDy*rho;
        flux    = -Dy*rho;        
        
        drhodt(bound) = rho(bound) - 1; 
        
        data.rho  = rho;
        data.flux = flux;
        drhodt    = MI.ApplyIntersectionBCs(drhodt,data,BCs);
        
    end

    function drhodt = rhsNoFlux(t,rho)
                                
        drhodt  = DDy*rho;
        flux    = -Dy*rho;        
        
        drhodt(bound) = flux(bound);
        
        data.rho  = rho;
        data.flux = flux;
        drhodt    = MI.ApplyIntersectionBCs(drhodt,data,BCs);
        
    end

    function rhoE = ExactSolutionPrescribed(y,t)
        rhoE = sqrt(t0/t)*exp(-y.^2/(4*t));
    end

    function rhoE = ExactSolutionConstant(y,t)
        z = (y-yMin)/(yMax-yMin)*2*pi;
        alpha = (2*pi/(yMax-yMin))^2;
        rhoE = sin(z) * exp(-alpha*t) + 1;
        % plus one for relative error to be better behaved
    end

    function rhoE = ExactSolutionNoFlux(y,t)
        z = (y-yMin)/(yMax-yMin)*2*pi;
        alpha = (2*pi/(yMax-yMin))^2;
        rhoE = cos(z) * exp(-alpha*t) + 1;
        % plus one for relative error to be better behaved
    end


end