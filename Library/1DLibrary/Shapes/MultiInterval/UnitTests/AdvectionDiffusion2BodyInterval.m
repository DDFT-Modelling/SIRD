function [line,Rho_t] = AdvectionDiffusion2BodyInterval(opts)

    if(nargin==0)
        geom.N = 10;
        geom.yMin = -1;
        geom.yMax = 1;
    else
        geom.N = opts.N;
        geom.yMin = opts.endpoints(1);
        geom.yMax = opts.endpoints(end);
    end

    line = SpectralLine(geom);
    
    
    [Pts,Diff,Int,Ind] = line.ComputeAll();    
    Interp             = line.ComputeInterpolationMatrix((-1:0.02:1)',true);                      
    
    Conv = line.ComputeConvolutionMatrix(@Gaussian,true);
    
    Dy = Diff.Dy;
    bound = Ind.bound;
    normal = Ind.normal;
    
    rho_ic = InitialCondition(Pts.y,0);
    
    D0 = 1;
    
    tMax = 0.5;
    
    outTimes = [0:tMax/100:tMax];
    
    mM        = ones(size(Pts.y));
    mM(bound) = 0;
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    [~,Rho_t] = ode15s(@rhs,outTimes,rho_ic,opts);

    maxRho = max(max(Rho_t));
    
    figure
    
    for iTimes = 1:length(outTimes)
        rhot = Rho_t(iTimes,:)';
        
        hold off
        line.plot(rhot);
        
        ylim([0,maxRho]);
        
        pause(0.05)
        
    end
    
    %----------------------------------------------------------------------
    
    function drhodt = rhs(t,rho)
                
        flux = getFlux(rho);
        
        drhodt = Dy*flux;

        flux = getFlux(rho);
        drhodt(bound) = normal*flux;
        
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
        
        y = exp(-y.^2);
       
    end

    function g = Gaussian(y)
        alpha = 2;
        g =  alpha * exp(- y.^2);
    end

end