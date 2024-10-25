 classdef FluctuatingChannel < M1SpectralSpectral
    properties 
        L1

        ft
        fb
        ftParams
        fbParams
        
    end
    
    methods
        function this = FluctuatingChannel(Geometry)
            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.L1= Geometry.L1;
            if(isfield(Geometry,'Origin'))              
                this.Origin = Geometry.Origin;
            end
            
            this.ft = str2func(Geometry.ft);
            this.fb = str2func(Geometry.fb);
            this.ftParams = Geometry.ftParams;
            this.fbParams = Geometry.fbParams;
            
            
            InitializationPts(this);            
            
            this.polar = 'cart';
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            L1      = this.L1;            
            ft      = this.ft;
            fb      = this.fb;
                  
            ftParams = this.ftParams;
            fbParams = this.fbParams;
            
            n  = length(x1);
            
            [y1_kv,dy1dx1] =  LinearMap(x1,0,L1);
            
            [top, dftdx1, dftdxx1] = ft(x1,ftParams); 
            [bottom, dfbdx1, dfbdxx1] = fb(x1,fbParams);

            
            [y2_kv,dy2dx2] =  LinearMap(x2,bottom,top);
            
            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = dy1dx1;
                %J(:,1,2) = zeros(n,1);
                J(:,2,2) = dy2dx2; % = dftdx1-dfbdx1/2
                J(:,2,1) = dfbdx1 + (dftdx1-dfbdx1).*(x2+1)/2;
            end

            if(nargout >= 4)
                % Hessian of y1 map
                dH1        = zeros(n,2,2);     
            end

            if(nargout >= 4)
                % Hessian of y2 map
                % y2 = bottom(x1) + (x2+1).*(top(x1)-bottom(x1))/2; 
                dH2        = zeros(n,2,2);            
                dH2(:,1,1) = dfbdxx1 + (dftdxx1-dfbdxx1).*(x2+1)/2;
                dH2(:,1,2) = (dftdx1-dfbdx1)/2;
                dH2(:,2,1) = (dftdx1-dfbdx1)/2;
            end

 

        end
        function [x1,x2] = CompSpace(this,y1,y2)
            exc = MException('FluctuatingChannel:CompSpace','not yet implemented');
            throw(exc);
        end    
        
        function Ind    = ComputeIndices(this)
            Ind      = GetIndicesFluctuatingChannel(this);
            this.Ind = Ind;
        end    
        
        
        
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@M1SpectralSpectral(this);
            %Check Accuracy            
%             if(this.sphere)
%                 y1s  = this.Pts.y1_kv;
%                 y2s  = this.Pts.y2_kv;
%                 int  = 2*int.*real(sqrt(this.R^2-y1s.^2-y2s.^2))';  
%                 this.Int = int;
% 
%                 ht   = this.R - this.h;
%                 area = pi*ht^2/3*(3*this.R - ht);
%             else
%                 th   = 2*acos(this.h/this.R);
%                 area = this.R^2/2*(th-sin(th));
%             end
%             if(nargout < 2)
%                 if(area == 0)
%                     disp(['Segment: Area is zero, Absolute error is: ',...
%                                     num2str(area-sum(this.Int))]);
%                 else
%                     disp(['Segment: Error of integration of area (ratio): ',...
%                                         num2str(1-sum(this.Int)/area)]);
%                 end
%             end
        end
            
    end
end