classdef Box < SpectralSpectral

    properties        
        y1Min = 0
        y2Min = 0
        y1Max,y2Max               
    end
    
    methods        
        function this = Box(Geometry)
            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            if(isfield(Geometry,'L1'))
                if(isfield(Geometry,'Origin'))
                    this.y1Min = Geometry.Origin(1);
                    this.y2Min = Geometry.Origin(2);
                end
                this.y1Max      = this.y1Min + Geometry.L1; 
                this.y2Max      = this.y2Min + Geometry.L2;            
            else
                this.y1Min = Geometry.y1Min;
                this.y1Max = Geometry.y1Max;                
                this.y2Min = Geometry.y2Min;
                this.y2Max = Geometry.y2Max;
            end
            
            this.polar = 'cart';
            
            InitializationPts(this);            
        end                         
    end
    
    methods (Access = public)
        
        function M_conv = ComputeConvolutionMatrix_Pointwise(this,f,shapeParams)
            % Convolution \int_D(x) f(x-y)g(y) dy
            
            disp('Computing Convolution matrices...'); 
                            
            
            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            if(isfield(shapeParams,'R'))                
                geom.N = shapeParams.N;
                geom.R = shapeParams.R;
            else
                % whole box
                M_conv = ComputeConvolutionMatrix(this,f);
                return;
                
            end

            if(isfield(shapeParams,'NT'))
                NT = shapeParams.NT;
            else
                NT = [10;10];
            end
                
            Pts = this.Pts;
            
            if(useDistance)
                fPTemp = f(GetDistance(this,Pts.y1_kv,Pts.y2_kv));
            else
                fPTemp = f(Pts.y1_kv,Pts.y2_kv);
            end
            
            fDim = size(fPTemp,2);
            
            N1  = this.N1;  N2  = this.N2;
            N1N2 = N1*N2;
            M_conv = zeros(N1*N2,N1*N2,fDim);
            
            wb = waitbar(0,'Computing Convolution Matrices');
            
            y1 = Pts.y1_kv;
            y2 = Pts.y2_kv;
            
            for i=1:N1N2
                
                waitbar(i/N1N2,wb);

                geom.Origin = [y1(i); y2(i)];
                disc = Disc(geom);
                disc.NT = NT;
                subShape = Intersect(this,disc);
                
                subConvShapePts = subShape.pts;
                
                % the interpolation onto D(x) from the full domain
                IP_i     = SubShapePts(this,subConvShapePts);
                
                % the integration weights for the subshape
                Int_i    = subShape.int;
                
                % f(x-y)
                dy1 = y1(i) - subShape.pts.y1_kv;
                dy2 = y2(i) - subShape.pts.y2_kv;
                
                if(useDistance)
                    fP          = f(GetDistance(this,dy1,dy2));
                else
                    fP          = f(dy1,dy2);
                end
                
                for iF = 1:fDim
                    M_conv(i,:,iF) = (Int_i.*fP(:,iF)')*IP_i;
                end

            end

            M_conv(isnan(M_conv)) = 0;
                
            disp('...done.');
            
            close(wb);
            
        end % convolution               

               
        
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)            
            [y1,dy1,dx,ddx,dddx,ddddx] = LinearMap(x1,this.y1Min,this.y1Max);
        end
        function xf = CompSpace1(this,y1)
            xf  = InvLinearMap(y1,this.y1Min,this.y1Max);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.y2Min,this.y2Max);
        end
        function x2 = CompSpace2(this,y2)                        
            x2  = InvLinearMap(y2,this.y2Min,this.y2Max);
        end                        
        
    end
end