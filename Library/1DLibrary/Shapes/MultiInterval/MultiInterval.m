classdef MultiInterval < handle
    
    properties (Access = public)                

        Intervals
        
        Pts
        Diff
        Ind
        Int
        Conv
        
        Interp
        
        Intersections
        nIntervals
        M

        yMask
        
    end
   
    methods
        
        function this = MultiInterval(intervals)

            nIntervals = length(intervals);
            this.nIntervals = nIntervals;

            for iInterval = 1:nIntervals
                intervalFn = str2func(intervals(iInterval).interval);
                this.Intervals(iInterval).Interval = intervalFn(intervals(iInterval).geom);
                this.Intervals(iInterval).Geom = intervals(iInterval).geom;
            end
            
            GetPoints(this);
            GetShapeMasks(this);

            ComputeDifferentiationMatrices(this);
            ComputeIndices(this);
            ComputeIntegrationVectors(this);
            ComputeInterpolationMatrices(this);
        end
        
        
        function GetPoints(this)
            ptsFields = {'x','y'};
            for iField = 1:length(ptsFields)
                allPts = [];
                for iInterval = 1:this.nIntervals
                    intervalPts = this.Intervals(iInterval).Interval.Pts.(ptsFields{iField});
                    allPts = [allPts; intervalPts]; %#ok
                end
                this.Pts.(ptsFields{iField}) = allPts;
            end
                
            this.Pts.yMin = min(min(this.Pts.y));
            this.Pts.yMax = max(max(this.Pts.y));

        end
        
        function GetShapeMasks(this)
            
            N = zeros(this.nIntervals,1);
            
            for iInterval = 1:this.nIntervals
                N(iInterval) = this.Intervals(iInterval).Interval.Pts.N;
            end
            
            this.M = sum(N);
            
            ptsCount = 0;
            
            for iInterval = 1:this.nIntervals
                mask = false(this.M,1);
                mask(ptsCount + 1 : ptsCount + N(iInterval)) = true;
                this.Intervals(iInterval).PtsMask = mask;
                ptsCount = ptsCount + N(iInterval);                
            end
            
        end
        
        
        function ComputeDifferentiationMatrices(this)

            for iInterval = 1:this.nIntervals
                this.Intervals(iInterval).Interval.ComputeDifferentiationMatrix;
            end
            
            diffFields = {'Dy','DDy'};
            
            for iField = 1:length(diffFields)
                allD = [];
                for iInterval = 1:this.nIntervals
                    intervalD = this.Intervals(iInterval).Interval.Diff.(diffFields{iField});
                    allD = blkdiag(allD, intervalD);
                end
                this.Diff.(diffFields{iField}) = allD;
            end
                            
        end
        
        function ComputeInterpolationMatrices(this,interpGrid)
            
            if(nargin<2)
                interpGrid = (-1:0.01:1)';
            end
                 
            InterPol = [];
            pts = [];
            
            for iInterval = 1:this.nIntervals
                IPi = this.Intervals(iInterval).Interval.ComputeInterpolationMatrix(interpGrid,true); 
                InterPoli = IPi.InterPol;
                ptsi = IPi.pts;
                
                InterPol = blkdiag(InterPol,InterPoli);
                pts = [pts;ptsi]; %#ok
            end
            
            this.Interp.InterPol = InterPol;
            this.Interp.pts = pts;
            
        end
        
        function IPU = ComputeInterpolationMatricesUniform(this,dy)
            
            if(~isfield(this.Pts,'Uniform'))
                if(nargin<2)
                    this.ComputeUniformPts;
                else
                    this.ComputeUniformPts(dy);
                end
            end
    
            
            IPU = [];
            
            for iInterval = 1:this.nIntervals
                if(nargin<2)
                    IPUi = this.Intervals(iInterval).Interval.InterpolationPlotUniform();
                else
                    IPUi = this.Intervals(iInterval).Interval.InterpolationPlotUniform(dy);
                end
                
                IPU = blkdiag(IPU,IPUi);
                
            end
            
            this.Interp.InterPolUniform = IPU;
            
        end
            
        function ComputeIndices(this)
            
            for iInterval = 1:this.nIntervals
                this.Intervals(iInterval).Interval.ComputeIndices;
            end
                                                
            % Compute intersections of left and right for each of the
            % intervals
            
            lrtbFields = {'left','right'};
            nFields = length(lrtbFields);
            
            for iInterval = 1:this.nIntervals
                                
                iInd = this.Intervals(iInterval).Interval.Ind;
                
                yi = this.Intervals(iInterval).Interval.Pts.y;
                
                for jInterval = 1:this.nIntervals
                    
                    if(iInterval ~= jInterval)

                        maskij = false(this.M,1);

                        jInd = this.Intervals(jInterval).Interval.Ind;

                        yj = this.Intervals(jInterval).Interval.Pts.y;

                        for iField = 1:nFields

                            mask1 = iInd.(lrtbFields{iField});
                            y1 = yi(mask1);

                            for jField = 1:nFields

                                mask2 = jInd.(lrtbFields{jField});
                                y2 = yj(mask2);              
                                
                                matches = this.CheckIntersection(y1,y2);

                                if(matches)
                                    
                                    maskij(this.Intervals(iInterval).PtsMask) = mask1;

                                    this.Intersections(iInterval,jInterval).PtsMask = maskij;

                                    this.Intersections(iInterval,jInterval).Side = lrtbFields{iField};

                                end

                            end
                                                        
                        end
                        
                    else % iInterval = jInterval
                    
                        this.Intersections(iInterval,jInterval).PtsMask = false(this.M,1);
                        this.Intersections(iInterval,jInterval).Side = [];
                    end
 
                    % Not necessary for 1D, but allows reusing some 2D code
                    this.Intersections(iInterval,jInterval).Flip = false;
                    
                end
                
            end
            
            % Find boundaries of MultiInterval by taking the union of all of
            % the corresponding boundaries
            
            boundaries = [];

            for iInterval = 1:this.nIntervals

                % Start with all edges being on the boundary of the
                % MultiInterval
                boundariesTemp = lrtbFields;
                for jInterval = 1:this.nIntervals
                    
                    % Find edges that aren't intersections with other
                    % intervals
                    intersectSide = this.Intersections(iInterval,jInterval).Side;
                    notIntersect = ~strcmp(intersectSide,boundariesTemp);
                    
                    % Keep only these edges
                    boundariesTemp = boundariesTemp(notIntersect);

                end
                
                % Construct contribution to the boundary from
                % Interval(iInterval)
                boundaries = cat(2,boundaries,{boundariesTemp});

            end
            
            % Loop through left, right combine these
            % edges that come from each of Interval(iInterval) into the
            % correseponding edge for the MultiInterval
            for iField = 1:nFields
               
                boundary = false(this.M,1);
                
                % the current edge
                field = lrtbFields{iField};
                
                for iInterval = 1:this.nIntervals
                    
                    % If the edge of Interval(iInterval) is in the boundary then
                    % add that to the correseponding boundary of MultiInterval
                    if(any(strcmp(field,boundaries{iInterval})))
                        
                        % need to ensure that the edge mask goes into the
                        % corresponding part of the MultiInterval pts by using
                        % the PtsMask for Intervals(iInterval)
                        boundary(this.Intervals(iInterval).PtsMask) = this.Intervals(iInterval).Interval.Ind.(field);
                        
                    end
                    
                end
                
                % Save mask for MultiInterval
                this.Ind.(field) = boundary;
            end
            
            % Compute entire boundary
            this.Ind.bound = this.Ind.left | this.Ind.right;

            % Easy indexing of subshape boundaries
            this.Ind.Shape = struct;
            for iField = 1:nFields
               
                % the current edge
                field = lrtbFields{iField};
                
                for iInterval = 1:this.nIntervals
                    % need to ensure that the edge mask goes into the
                    % corresponding part of the MultiInterval pts by using
                    % the PtsMask for Interval(iInterval)
                    edge = false(this.M,1);
                    edge(this.Intervals(iInterval).PtsMask) = this.Intervals(iInterval).Interval.Ind.(field);
                    this.Ind.Interval(iInterval).(field) = edge;
                end
                                
            end
                                    
            % Construct intersection mask
            intersections = false(this.M,1);
            for iInterval = 1:this.nIntervals
                for jInterval = 1:this.nIntervals
                    intersections(this.Intersections(iInterval,jInterval).PtsMask) = true;
                end
            end
            this.Ind.intersections = intersections;
                            
            % Construct normals for the MultiInterval
            boundFields = {'left','right','bound'};
            normalFields = {'normalLeft','normalRight','normal'};
            nFields = length(boundFields);
            
            for iField = 1:nFields

                allNormals = [];

                for iInterval = 1:this.nIntervals

                   Ni = this.Intervals(iInterval).Interval.Pts.N;

                   bound = this.Intervals(iInterval).Interval.Ind.(boundFields{iField}); 

                   normal = this.Intervals(iInterval).Interval.Ind.(normalFields{iField}); 

                   temp = zeros(Ni,Ni);
                   temp(bound,:) = normal;

                   allNormals = blkdiag(allNormals,temp);

                end

                allNormals(~this.Ind.(boundFields{iField}),:) = []; %#ok
                this.Ind.(normalFields{iField}) = sparse(allNormals);
            
            end
           
            % Construct normals for intersections
            for iField = 1:nFields
                
                for iInterval = 1:this.nIntervals
                    for jInterval = 1:this.nIntervals
                        
                        if(iInterval ~= jInterval)
                        
                            sidei = this.Intersections(iInterval,jInterval).Side;
                            iPos = find(strcmp(sidei,boundFields));

                            if(~isempty(iPos))

                                temp = zeros(1,this.M);
                                %maski = this.Intersections(iInterval,jInterval).PtsMask;
                                maski = this.Intervals(iInterval).PtsMask;
                                
                                normali = this.Intervals(iInterval).Interval.Ind.(normalFields{iPos});
                                temp(1,maski) = normali;
                                this.Intersections(iInterval,jInterval).Normal = sparse(temp);
                            end

                        end
                        
                    end
                end
            end
            
        end
        
        function dydt = ApplyIntersectionBCs(this,dydt,data,BCs)
        
            for iInterval = 1:this.nIntervals
                for jInterval = iInterval+1:this.nIntervals

                    if(~isempty(BCs(iInterval,jInterval).function))
                        BCfn = str2func(BCs(iInterval,jInterval).function);
                        [BC,mask] = BCfn(data,this.Intersections,iInterval,jInterval);
                        dydt(mask) = BC(mask);
                    end                
                end
            end            
            
        end
        
        function ComputeIntegrationVectors(this)
            
            this.Int = [];
            
            for iInterval = 1:this.nIntervals
                Inti = this.Intervals(iInterval).Interval.ComputeIntegrationVector;
                this.Int = [this.Int, Inti];
            end
        end
                
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            y = this.Pts.y;
            
            fPTemp = f(y);
            
            fDim = size(fPTemp);
            
            nElts = prod(fDim(2:end));
            
            IntT = (this.Int).';  % M x 1
            
            IntT = IntT(:,ones(1,nElts)); % M x nElts
            IntT = reshape(IntT,fDim);    % size(f)
            
            M_conv = zeros([this.M,this.M,fDim(2:end)]);
            
            Mmask = repmat({':'},[1,fDim]);
            
            for i=1:this.M 
                fP          = f(y(i) - y);
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            M_conv(isnan(M_conv)) = 0;
            
            if((nargin >= 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
            
        end                
        
        %------------------------------------------------------------------
        % Plotting
        %------------------------------------------------------------------
        
        function PlotGrid(this)
            for iInterval = 1:this.nIntervals
                thisInterval = this.Intervals(iInterval).Interval;
                thisInterval.PlotGrid;
                hold on
            end 
        end
    
        function PlotBound(this,boundList)
            if(nargin==1)
                boundList = {'bound'};
            end
            
            colourList = {'r','b','m','c'};
            
            y = this.Pts.y;
            z = zeros(size(y));
            
            nFields = length(boundList);
            
            for iField = 1:nFields
                mask = this.Ind.(boundList{iField});
                scatter(y(mask),z(mask),'o','MarkerEdgeColor',colourList{iField},'MarkerFaceColor',colourList{iField});
                hold on
            end
            legend(boundList);
        end
        
        function PlotIntersections(this,plotIJ,opts)
            
            y = this.Pts.y;
            z = zeros(size(y));
            
            if(nargin<2 || isempty(plotIJ))
                plotIJ = true(this.nIntervals);
            end
            
            if(nargin<3)
                opts = [];
            end
            
            if(isfield(opts,'colour'))
                colour = opts.colour;
            else
                colour = 'b';
            end
                        
            for iInterval = 1:this.nIntervals
                for jInterval = iInterval+1:this.nIntervals
                    if(plotIJ(iInterval,jInterval))
                        PtsMask = this.Intersections(iInterval,jInterval).PtsMask;            
                        scatter(y(PtsMask),z(PtsMask),'x',colour);
                        hold on
                    end
                end
            end
            
            hold off
              
            this.SetAxes;
                       
        end
            
        function PlotNormalsIntervals(this)
            
            for iInterval = 1:this.nIntervals
                y = this.Intervals(iInterval).Interval.Pts.y;
                z = zeros(size(y));
                
                dir1 = ones(size(y));
                
                normal = this.Intervals(iInterval).Interval.Ind.normal;
                bound = this.Intervals(iInterval).Interval.Ind.bound;
                
                v1 = normal * dir1;
                v2 = zeros(size(v1));
                
                vnorm = sqrt(v1.^2 + v2.^2);
                v1 = v1./vnorm;
                v2 = v2./vnorm;
                
                y1 = y(bound);
                y2 = z(bound);
                
                quiver(y1, y2, v1, v2, 'AutoScale','off');
                
                hold on
                
            end
       
            axis equal
            
        end

        
        function PlotNormals(this,boundList)
            
            boundListFull = {'left','right','bound'};
            normalListFull = {'normalLeft','normalRight','normal'};
            
            if(nargin==1)
                boundList = {'bound'};
            end
            
            colourList = {'r','b','m','c','k'};

            y = this.Pts.y;
            z = zeros(size(y));
            
            dir1 = ones(size(y));
            
            nFields = length(boundList);
            
            for iField = 1:nFields

                iPos = find(strcmp(boundList{iField},boundListFull));
                
                normal = this.Ind.(normalListFull{iPos});
                bound = this.Ind.(boundListFull{iPos});

                v1 = normal * dir1;
                v2 = zeros(size(v1));

                vnorm = sqrt(v1.^2 + v2.^2);
                v1 = v1./vnorm;
                v2 = v2./vnorm;

                y1Plot = y(bound);
                y2Plot = z(bound);         

                quiver(y1Plot, y2Plot, v1, v2, 'AutoScale','off','Color',colourList{iPos});                
                
                hold on
            end

            axis equal
            
        end
     
        function PlotNormalsIntersections(this)
            
            y = this.Pts.y;
            z = zeros(size(y));
            
            dir1 = ones(size(y));
            
            for iInterval = 1:this.nIntervals
                for jInterval = 1:this.nIntervals
               
                    if(any(this.Intersections(iInterval,jInterval).PtsMask))
                        normal = this.Intersections(iInterval,jInterval).Normal;
                        bound = this.Intersections(iInterval,jInterval).PtsMask;

                        v1 = normal * dir1;
                        v2 = zeros(size(v1));

                        vnorm = sqrt(v1.^2 + v2.^2);
                        v1 = v1./vnorm;
                        v2 = v2./vnorm;
                        
                        y1Plot = y(bound);
                        y2Plot = z(bound);         

                        quiver(y1Plot, y2Plot, v1, v2, 'AutoScale','off');

                        hold on
                    end
                    
                end
                
            end
            
            axis equal
            
        end
        
        function Plot(this,V,opts)

            if(nargin<3)
                opts.linecolor = 'k';
                opts.linewidth = 2;
            end
            
            
            for iInterval = 1:this.nIntervals
                Vi = V(this.Intervals(iInterval).PtsMask);
                this.Intervals(iInterval).Interval.plot(Vi,opts);
                hold on
            end
            
            this.SetAxes;
            
        end
        

        function SetAxes(this,ha)
            
            if (nargin < 2)
                ha = gca;
            end
            
            xMin = this.Pts.yMin;
            xMax = this.Pts.yMax;
            %axes(ha);
            xlim([xMin, xMax]);
            
            %pbaspect(ha,[(xMax-xMin) (yMax-yMin) 1/2*min((xMax-xMin),(yMax-yMin))]);
        end

        function UniformPts = ComputeUniformPts(this,dy)
            
            y = [];
            
            for iInterval = 1:this.nIntervals
                if(nargin<2)
                    Ptsi = this.Intervals(iInterval).Interval.ComputeUniformPts;
                else
                    Ptsi = this.Intervals(iInterval).Interval.ComputeUniformPts(dy);
                end
                
                y = [y; Ptsi.y]; %#ok
            end
                
            % remove points from intersections which may be double counted
            %UniformPts.y = unique(y); 
            UniformPts.y = y; 
            
            this.Pts.Uniform = UniformPts;
            
        end
        
        function matches = myisequal(this,y1,y2)
            
            matches = 0;
            
            if(size(y1) == size(y2))
                d = abs(y1-y2);
                if max(d)< 10^(-14)
                    matches = 1;
                end
            end
        end
        
        function matches = CheckIntersection(this,y1,y2)            
            matches = this.myisequal(y1,y2);
        end
        
    end % end methods
    
end