function MI = MI_TanhTest()

    close all

    if(nargin==0)
        opts = struct;
    end
        
    if(isfield(opts,'endpoints'))
        endpoints = opts.endpoints;
    else
        %endpoints = [-1;-0.5;0.5;1];
        endpoints = [-1;1];
    end

    nIntervals = length(endpoints)-1;
    
    if(isfield(opts,'N'))
        N = opts.N;
    else
        %N = [10;100;10];
        N = 100;
    end

    if(length(N)==1)
        N = repmat(N,nIntervals,1);
    end
    
    if(isfield(opts,'alpha1'))
        alpha1 = opts.alpha1;
    else
        alpha1 = 10;
    end
    
    if(length(alpha1)==1)
        alpha1 = repmat(alpha1,nIntervals,1);
    end

    if(isfield(opts,'alpha1'))
        alpha2 = opts.alpha1;
    else
        alpha2 = 0;
    end
    
    if(length(alpha2)==1)
        alpha2 = repmat(alpha2,nIntervals,1);
    end
    
    if(isfield(opts,'interval'))
        interval = opts.interval;
    else
        %interval = {'SpectralLine','SpectralLineTan','SpectralLine'};
        interval = {'SpectralLineTan'};
    end
    
    if(isfield(opts,'doPlots'))
        doPlots = opts.doPlots;
    else
        doPlots = true;
    end
    
    intervals(nIntervals).interval = []; % preallocate
    
    for iInterval = 1:nIntervals
        intervals(iInterval).interval = interval{iInterval};
        geom.N    = N(iInterval);
        geom.yMin = endpoints(iInterval);
        geom.yMax = endpoints(iInterval+1);
        geom.alpha1 = alpha1(iInterval);
        geom.alpha2 = alpha2(iInterval);
        intervals(iInterval).geom = geom;
    end
    
    MI = MultiInterval(intervals);
    y = MI.Pts.y;
    Diff = MI.Diff;

    figure('Position',[100,100,1200,800]);
    %subplot(2,1,1)
    %MI.PlotGrid;

    alpha = 50;
    beta = 1;

    [dfe,ddfe] = diffTest(y,beta,Diff);

    disp(['Errors x^4: D1=' num2str(dfe) ', D2=' num2str(ddfe)]);

    [dfe,ddfe] = diffTest2(y,alpha,Diff);

    disp(['Errors tanh: D1=' num2str(dfe) ', D2=' num2str(ddfe)]);

    f4 = y4(y,beta);

    subplot(3,1,1)
    MI.Plot(f4);

    ft = t(y,alpha);

    hold on
    MI.Plot(ft);
    
    subplot(3,1,2)
    MI.Plot(Diff.Dy*ft);
    subplot(3,1,3)
    MI.Plot(Diff.DDy*ft);

    function [f,df,ddf] = y4(y,beta)
        z = y/beta;
        f     = z.^4;
        df    = 4*z.^3/beta;
        ddf   = 12*z.^2/beta^2;
    end

    function [f,df,ddf] = y4D(y,beta,Diff)
        f     = (y/beta).^4;
        df    = Diff.Dy*f;
        ddf   = Diff.DDy*f;
    end

    function [dfe,ddfe] = diffTest(y,beta,Diff)
        [f,df,ddf] = y4(y,beta);
        [f,Df,DDf] = y4D(y,beta,Diff);
        dfe = max(abs(df-Df)./abs(df+1));
        ddfe = max(abs(ddf-DDf)./abs(ddf+1));
    end

    function [f,df,ddf] = t(y,alpha)
        f     = tanh(alpha*y);
        df    = alpha*sech(alpha*y).^2;
        ddf    = -2*alpha^2*tanh(alpha*y).*sech(alpha*y).^2;
    end

    function [f,df,ddf] = tD(y,alpha,Diff)
        f     = tanh(alpha*y);
        df    = Diff.Dy*f;
        ddf   = Diff.DDy*f;
    end

    function [dfe,ddfe] = diffTest2(y,alpha,Diff)
        [f,df,ddf] = t(y,alpha);
        [f,Df,DDf] = tD(y,alpha,Diff);
        dfe = max(abs(df-Df)./abs(df+1));
        ddfe = max(abs(ddf-DDf)./abs(ddf+1));
    end


end