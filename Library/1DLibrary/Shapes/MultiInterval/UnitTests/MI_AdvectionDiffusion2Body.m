function [MI,Rho_t] = MI_AdvectionDiffusion2Body(opts)

    if(nargin==0)
        opts = struct;
    end
        
    if(isfield(opts,'endpoints'))
        endpoints = opts.endpoints;
    else
        endpoints = [-1;-0.5;0.5;1];
        %endpoints = [-1;1];
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

    if(isfield(opts,'doPlots'))
        doPlots = opts.doPlots;
    else
        doPlots = true;
    end
    
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
    Int = MI.Int;
    Ind = MI.Ind;
    
    y = Pts.y;
    bound = Ind.bound;
    normal = Ind.normal;
    intersections = Ind.intersections;
    
    Conv = MI.ComputeConvolutionMatrix(@Gaussian,true);
    

    
    rho_ic = InitialCondition(y,0);
    
    D0 = 1;
    
    tMax = 0.5;
    
    outTimes = [0:tMax/100:tMax];
    
    mM        = ones(size(Pts.y));
    mM(bound) = 0;
    mM(intersections) = 0;
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    [~,Rho_t] = ode15s(@rhs,outTimes,rho_ic,opts);

    if(doPlots)
        maxRho = max(max(Rho_t));

        figure

        for iTimes = 1:length(outTimes)
            rhot = Rho_t(iTimes,:)';

            hold off
            MI.Plot(rhot);

            ylim([0,maxRho]);

            pause(0.05)

        end
    end
        
    %----------------------------------------------------------------------
    
    function drhodt = rhs(t,rho)
                
        flux = getFlux(rho);
        
        drhodt = Dy*flux;

        flux = getFlux(rho);
        drhodt(bound) = normal*flux;
        
        data.rho  = rho;
        data.flux = flux;
        drhodt    = MI.ApplyIntersectionBCs(drhodt,data,BCs);

        drhodt = drhodt(:);
                
    end

    function f = getFlux(rho)
        j = getAdvection(rho);
        f2 = get2Body(rho);
        f = D0*Dy*rho + j + f2;
    end

    function j = getAdvection(rho)
        v = ones(size(rho));
        j = - rho.*v;
    end

    function f = get2Body(rho)
        f = rho.*(Dy*(Conv * rho));
    end


    function y = InitialCondition(y,t)
        
        y = exp(-10*y.^2);
       
    end

    function g = Gaussian(y)
        alpha = 5;
        g =  alpha * exp(- y.^2);
    end

end