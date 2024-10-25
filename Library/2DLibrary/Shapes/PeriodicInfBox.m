classdef PeriodicInfBox < SpectralFourier

    properties       
        % y1 spectral inf, y2 periodic
        %y1Min = 0
        y2Min = 0
        %y1Max
        y2Max  
        y1_0 = 0; % this origin isn't implemented yet
        L1 % note that this isn't the size of the y1 domain; that's infinite, 
           % rather it's a measure of how spread out the points are
        L2
    end
    
    methods        
        function this = PeriodicInfBox(Geometry)
           this@SpectralFourier(Geometry.N(1),Geometry.N(2));
                        
            this.y2Min = Geometry.y2Min;
            this.y2Max = Geometry.y2Max;
            
            this.L1 = Geometry.L1;
            this.L2 = this.y2Max - this.y2Min;
            
            InitializationPts(this);      
            this.polar = 'cart';
        end                         
    
        % overridden from SpectralFourier to deal with infinite entries in
        % the integration vector
        function int = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@SpectralFourier(this);
            int(~isfinite(int)) = 0; % deal with points at \pm \infty
        end
        
        
        % overridden from SpectralFourier
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)

            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            N1  = this.N1;  
            N2  = this.N2;
            Pts = this.Pts;
            Int = this.Int;  % 1 x N1*N2

            if(useDistance)
                fPTemp = f(GetDistance(this,Pts.y1_kv,Pts.y2_kv));
            else
                fPTemp = f(Pts.y1_kv,Pts.y2_kv);
            end

            fDim = size(fPTemp);

            nElts = prod(fDim(2:end));

            IntT = Int.';  % N1*N2 x 1

            IntT = IntT(:,ones(1,nElts)); % N1*N2 x nElts
            IntT = reshape(IntT,fDim);    % size(f)

            M_conv = zeros([N1*N2,N1*N2,fDim(2:end)]);

            Mmask = repmat({':'},[1,fDim]);

            for i=1:(N1*N2) 
                if(useDistance)
                    fP          = f(GetDistance(this,Pts.y1_kv(i) - Pts.y1_kv,Pts.y2_kv(i) - Pts.y2_kv));
                else
                    fP          = f(GetDistance1(this,Pts.y1_kv(i) - Pts.y1_kv),GetDistance2(this,Pts.y2_kv(i) - Pts.y2_kv));
                    %fP          = f(Pts.y1_kv(i) - Pts.y1_kv,Pts.y2_kv(i) - Pts.y2_kv);
                end
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            M_conv(isnan(M_conv)) = 0;

            if((nargin >= 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
        end        


        % overridden from Shape
        function d  = GetDistance(this,pts_y1,pts_y2)             
            if(nargin == 2)
                pts_y2 = pts_y1.y2_kv;
                pts_y1 = pts_y1.y1_kv;
            end

            y1 = this.GetDistance1(pts_y1);   % spectral
            y2 = this.GetDistance2(pts_y2); % Fourier

            d       = sqrt(y1.^2 + y2.^2);            
        end

        function M_IntConv = ComputeIntegrationConvolutionMatrix(this,fInt,fConv,shapeParams)
            
            if(nargin<4)
                shapeParams = struct;
            end
            
            % Weighted integration in y1 (spectral)
            
            % set up spectral line
            geom1.N = this.N1;
%             geom1.yMin = this.y1Min;
%             geom1.yMax = this.y1Max;
            geom1.L = this.L1;
            sLine = InfSpectralLine(geom1);

            % compute integration weights
            Int1 = sLine.ComputeIntegrationVector;
            ws = sLine.Pts.y;

            gLine = fInt(ws);
            %TOM: removed non-finite points (11/08/20)
            gLine(~isfinite(gLine)) = 0;
            I1Line = Int1.* gLine';
            I1Rep = repmat(I1Line,[this.N1,1]);
            I1 = sparse(kron(I1Rep,eye(this.N2)));
            
            % Convolution in y2 (Fourier)
            
            % set up Fourier line
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            fLine = FourierLine(geom2);

            % compute convolution (which requires the integration vector)
            fLine.ComputeIntegrationVector;
            % format: C2 = fLine.ComputeConvolutionMatrix(this,f,shapeParams,parent)
            C2Line = fLine.ComputeConvolutionMatrix(fConv,shapeParams);
            
            C2 = sparse(kron(eye(this.N1),C2Line));
            
            % Full integral
            M_IntConv = C2*I1;
        end

        function [Int1,Int2] = Compute1DIntegrationMatrices(this)

            % Spectral integration
            geom1.N = this.N1;
%           TOM 11/08/2020
%             geom1.yMin = this.y1Min;
%             geom1.yMax = this.y1Max;
            geom1.L = this.L1;
            sLine = InfSpectralLine(geom1);
            I1 = sLine.ComputeIntegrationVector;
            Int1 = sparse(kron(I1,eye(this.N2)));
            
            % Fourier integration
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            fLine = FourierLine(geom2);
            I2 = fLine.ComputeIntegrationVector;
            Int2 = sparse(kron(eye(this.N1),I2));
        end
        
        function d  = GetDistance1(this,y1)             
            d = y1;
        end
        function d  = GetDistance2(this,y2)             
            d = mod(y2 + this.L2/2,this.L2) - this.L2/2;
        end

    end
    
    methods (Access = public)
        
        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            [y,dy,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);
        end
        function xf = CompSpace1(this,y1)
            xf  = InvSqrtMap(y1-this.y1_0,this.L1,inf);
        end               
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(x2,this.y2Min,this.y2Max);
        end
        function x2 = CompSpace2(this,y2)                        
            x2 = (y2-this.y2Min)/(this.y2Max-this.y2Min);   
        end        
        
    end
end