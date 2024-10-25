classdef Quadrilateral < M1SpectralSpectral
    properties
        Y % corner points, mapped in order from (-1,-1), (1,-1), (1,1), (-1,1)
          % Y(i,:) = (y_i1,y_i2)
        A % mapping matrix from (0,0), (1,0), (1,1), (0,1) -> Y
        
    end
    
    methods 
        function this = Quadrilateral(Geometry)
            
            Geometry.Y = SortQuadrilateralIndices(Geometry.Y);
            
            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.Y   = Geometry.Y;
            this.A   = GetQuadrilateralParameters(this.Y);
            this.polar = 'cart';
            InitializationPts(this);
        end
                
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            % see e.g. https://www.particleincell.com/2012/quad-interpolation/
            
            n  = length(x1);
            O  = ones(n,1);

            % shift [-1,1] to [0,1] via linear map
            x1 = (x1+1)/2;
            x2 = (x2+1)/2;
            
            Dx1 = 1/2;
            Dx2 = 1/2;
            
            % apply mapping from [0,1] x [0,1] to quadrilateral determined
            % by corners Y
            h = this.A*[O';x1';x2';(x1.*x2)'];
            y1_kv = h(1,:)';  y2_kv = h(2,:)';

            if(nargout >= 3)
                % Jacobian
                J        = zeros(n,2,2);
                J(:,1,1) = (this.A(1,2)*O + this.A(1,4)*x2) * Dx1;
                J(:,1,2) = (this.A(1,3)*O + this.A(1,4)*x1) * Dx2;
                J(:,2,1) = (this.A(2,2)*O + this.A(2,4)*x2) * Dx1;
                J(:,2,2) = (this.A(2,3)*O + this.A(2,4)*x1) * Dx2;
            end

            if(nargout >= 4)
                % Hessian of y1 map
                dH1        = zeros(n,2,2);     
                dH1(:,1,2) = this.A(1,4)*O * Dx1 * Dx2;
                dH1(:,2,1) = this.A(1,4)*O * Dx1 * Dx2;
            end

            if(nargout >= 4)
                % Hessian of y2 map
                dH2        = zeros(n,2,2);            
                dH2(:,1,2) = this.A(2,4)*O * Dx1 * Dx2;
                dH2(:,2,1) = this.A(2,4)*O * Dx1 * Dx2;
            end

        end
        function [x1,x2] = CompSpace(this,y1,y2)
            
            % map onto [0,1] x [0,1]
            alpha = this.A(1,:);
            beta = this.A(2,:);
            
            % a x2^2 + b x2 + c = 0
            a = alpha(4)*beta(3) - alpha(3)*beta(4);
            b = alpha(4)*beta(1) - alpha(1)*beta(4) + alpha(2)*beta(3) ...
                - alpha(3)*beta(2) + y1*beta(4) - y2*alpha(4);
            c = alpha(2)*beta(1) - alpha(1)*beta(2) + y1*beta(2) - y2*alpha(2);
            
            if(a==0)
                w2 = -c./b;
            else
                det = sqrt(b.^2-4*a.*c);
                w2 = (-b + det)/(2*a);
            end
            w1 = (y1 - alpha(1) - alpha(3)*w2)./(alpha(2) + alpha(4)*w2);
            
            % map onto [-1,1] x [-1,1]
            x1 = 2*w1 - 1;
            x2 = 2*w2 - 1;
        end
        
        
        function [int,area] = ComputeIntegrationVector(this)
            int  = ComputeIntegrationVector@M1SpectralSpectral(this);
            %Check Accuracy
            x = this.Y(:,1);
            y = this.Y(:,2);
            area = 1/2*abs( ( x(1)*y(2) - y(1)*x(2) ) + ( x(2)*y(3) - y(2)*x(3) ) ...
                             + ( x(3)*y(4) - y(3)*x(4) ) + ( x(4)*y(1) - y(4)*x(1) ));       
            
%              if(nargout < 2)
%                 if(area == 0)
%                     disp(['Quadrilateral: Error of integration of area(=0): ',...
%                                         num2str(area-sum(this.Int))]);
%                 else
%                     disp(['Quadrilateral: Error of integration of area(ratio): ',...
%                                         num2str(1-sum(this.Int)/area)]);
%                 end
%             end
        end
        
        function Ind    = ComputeIndices(this)
            Ind      = GetIndicesQuadrilateral(this);
            this.Ind = Ind;
        end    
                
    end
end

