function dataDisc = Intersect_Disc(MainShape,discShape)

    r   = discShape.R;
    N   = [discShape.N1,discShape.N2];
	y20 = discShape.Origin(2);       
    
	shape.N  = N;
    shape.R  = r;     
    shape.Origin = discShape.Origin;
    
    if((isa(discShape,'Disc') && (discShape.sphere == true)) || ...
             (isa(discShape,'Sphere') && discShape.volume == true))
            shape.volume = true;
    else
        shape.volume = false;
    end
    
   % fullDisc = false;
    
    if(isa(MainShape,'HalfSpace'))   
        
        y2Min = MainShape.y2Min;
        if(isa(MainShape,'HalfSpaceSkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end

        %1. find points of disc in HalfSpace            
        if(y20 == inf)
            shape.Origin(2) = 0;
            shape.theta1 = 0;
            shape.theta2 = pi;
           % fullDisc = true;
        elseif(y20 >= y2Min + r)
            %1a. if full disc is in HalfSpace                
            shape.theta1 = 0;
            shape.theta2 = pi;                             
        elseif((y20 < (y2Min + r)) && ...
               (y20 >= (y2Min - r)))
            %1b. if part of disc is in HalfSpace  (>= half)
            %1b1. Integrate over segment in HalfSpace
            %shape.Origin    = [0,y20];
            th                = acos((y20 - y2Min)/r);
            shape.theta1      = 0;
            shape.theta2      = pi-th;                                                
        %1c. if part of disc is in HalfSpace  (< half)    
        else
            exc = MException('HalfSpace_FMT:AverageDisc','case not implemented');
            throw(exc);                
        end            
        
    elseif(isa(MainShape,'Box'))
        area           = Intersect_Disc_Box(discShape,MainShape);                
        
    elseif(isa(MainShape,'InfCapillary') && (MainShape.y2Max - MainShape.y2Min)>= 2*r)   
%         
         y2Min = MainShape.y2Min;
         y2Max = MainShape.y2Max;
         if(isa(MainShape,'InfCapillarySkewed'))   
             y2Min = y2Min*sin(MainShape.alpha);
             y2Max = y2Max*sin(MainShape.alpha);
         end   
 
         %1. find points of disc fully in InfCapillary            
         if((y20 < (y2Min + r)) && (y20 >= (y2Min - r)))
             th                = acos((y20 - y2Min)/r);
             shape.theta1      = 0;
             shape.theta2      = pi-th;                                                
         elseif( (y20 >= y2Min + r) && (y20 <= y2Max - r))           
             shape.theta1 = 0;
             shape.theta2 = pi;                              
         elseif((y20 > (y2Max - r)) && (y20 <= (y2Max + r)))
             th                = acos((y2Max-y20)/r);
             shape.theta1      = th;
             shape.theta2      = pi;           
         else
             exc = MException('Intersect_Disc','case not implemented');
             throw(exc);                
         end            
    else
        exc = MException('Intersect_Disc','case not implemented');
        throw(exc);                
    end        
    
%     if(fullDisc)
%         shape.N  = N;
%         shape.R  = r;    
%         area = Disc(shape);
%     else
    if(~isa(MainShape,'Box'))
        area = Sphere(shape);
    end

    dataDisc.pts       = area.GetCartPts();           
    ptsLoc.y1_kv       = dataDisc.pts.y1_kv - shape.Origin(1);
    ptsLoc.y2_kv       = dataDisc.pts.y2_kv - shape.Origin(2);        
    dataDisc.ptsPolLoc = Cart2PolPts(ptsLoc);

    if(y20 == inf)
        dataDisc.pts.y2_kv = dataDisc.pts.y2_kv + y20;
    end

    [dataDisc.int,dataDisc.area]     = area.ComputeIntegrationVector();                                   
    
end   