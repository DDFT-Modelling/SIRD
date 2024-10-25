classdef BoxFourierFourier < FourierFourier
    properties        
        y1Min = 0
        y2Min = 0
        y1Max,y2Max 
        L1, L2
    end
    
    methods        
        function this = BoxFourierFourier(Geometry)
            this@FourierFourier(Geometry.N(1),Geometry.N(2));
            
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
            
            this.L1 = this.y1Max - this.y1Min;
            this.L2 = this.y2Max - this.y2Min;
            
            this.polar = 'cart';
            
            InitializationPts(this);            
        end                         
    end
    
    methods (Access = public)
        
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(x1,this.y1Min,this.y1Max);
        end
        function x1 = CompSpace1(this,y1)                        
            %x1 = (y1-this.y1Min)/(this.y1Max-this.y1Min);   
            x1 = InvLinearMap01(y1,this.y1Min,this.y1Max);
        end        

        
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(x2,this.y2Min,this.y2Max);
        end
        function x2 = CompSpace2(this,y2)                        
            %x2 = (y2-this.y2Min)/(this.y2Max-this.y2Min);   
            x2 = InvLinearMap01(y2,this.y2Min,this.y2Max);
        end        
                     
        
    end
end