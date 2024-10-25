classdef MultiShape < handle
    
    properties (Access = public)
        
        Shapes
        
        Pts
        Diff
        Ind
        Int
        Conv
        
        Interp
        
        Intersections
        nShapes
        M
        
        y1Mask
        y2Mask
        
    end
    
    methods
    
        %------------------------------------------------------------------
        % Initialisation
        %------------------------------------------------------------------

        function this = MultiShape(shapes,normalsOverride)
            
            % normalsOverride allows individual normals to be overriden,
            % for example when then intersection of shapes is not
            % perpendicular to the boundary
            if(nargin<2)
                normalsOverride = [];
            end
            
            nShapes = length(shapes);
            this.nShapes = nShapes;
            
            % construct shapes from input structure
            for iShape = 1:nShapes
                shapeFn = str2func(shapes(iShape).shape);
                this.Shapes(iShape).Shape = shapeFn(shapes(iShape).geom);
                this.Shapes(iShape).Geom = shapes(iShape).geom;
            end
            
            % construct standard quantities
            GetPoints(this);
            GetShapeMasks(this);
            ComputeDifferentiationMatrices(this);
            ComputeIndices(this,normalsOverride);
            ComputeIntegrationVectors(this);
            ComputeInterpolationMatrices(this);
        end

        %------------------------------------------------------------------
        % Physical and Computational points
        %------------------------------------------------------------------

        function GetPoints(this)
        % *** Note that all points in MS are in Cartesian coordinates ***

            % computational points xi_kv and physical points 
            ptsFields = {'x1_kv','x2_kv','y1_kv','y2_kv'};

            % collect points from each shape into a single vector
            for iField = 1:length(ptsFields)
                allPts = [];
                for iShape = 1:this.nShapes
                    shapePts = this.Shapes(iShape).Shape.Pts.(ptsFields{iField});
                    allPts = [allPts; shapePts]; %#ok
                end
                this.Pts.(ptsFields{iField}) = allPts;
            end
            
            % masks to allow easy access to coordinates on a Shapewise
            % basis
            this.GetShapeMasks();
            
            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;
            
            % convert from polar to Cartesian coordinates where necessary
            for iShape=1:this.nShapes
                if strcmp(this.Shapes(iShape).Shape.polar,'polar')
                    
                    % points are currently in polar
                    shapeMask = this.Shapes(iShape).PtsMask;
                    y1Polars = y1(shapeMask);
                    y2Polars = y2(shapeMask);

                    % convert to Cartesian using function from Shape
                    ptsCart = GetCartPts(this.Shapes(iShape).Shape,y1Polars,y2Polars);
                    y1Cart = ptsCart.y1_kv;
                    y2Cart = ptsCart.y2_kv;

                    % assign back to correct part of MS Pts
                    this.Pts.y1_kv(shapeMask) = y1Cart;
                    this.Pts.y2_kv(shapeMask) = y2Cart;
                    
                end
            end
            
            % determine points bounds
            this.Pts.y1Min = min(min(this.Pts.y1_kv));
            this.Pts.y1Max = max(max(this.Pts.y1_kv));
            this.Pts.y2Min = min(min(this.Pts.y2_kv));
            this.Pts.y2Max = max(max(this.Pts.y2_kv));
            
        end
        
        %------------------------------------------------------------------
        % Shapewise masks for scalars and vectors
        %------------------------------------------------------------------

        function GetShapeMasks(this)
            
            % number of points in each shape
            N1 = zeros(this.nShapes,1);
            N2 = zeros(this.nShapes,1);
            
            for iShape = 1:this.nShapes
                N1(iShape) = this.Shapes(iShape).Shape.Pts.N1;
                N2(iShape) = this.Shapes(iShape).Shape.Pts.N2;
            end
            
            Mi = N1.*N2;
            this.M = sum(Mi);
            
            % determine masks for each shape by iterating through

            % counters for how many points we've assigned so far
            % indices go shape-by-shape so simple iteration works here
            ptsCount = 0;
            ptsCountVec = 0;
            
            for iShape = 1:this.nShapes

                % scalar mask for each shape
                mask = false(this.M,1);
                mask(ptsCount + 1 : ptsCount + Mi(iShape)) = true;
                this.Shapes(iShape).PtsMask = mask;
                ptsCount = ptsCount + Mi(iShape);
                
                % vector mask for each shape
                maskVec = false(2*this.M,1);
                maskVec(ptsCountVec + 1 : ptsCountVec + 2*Mi(iShape)) = true;
                this.Shapes(iShape).PtsMaskVec = maskVec;
                ptsCountVec = ptsCountVec + 2*Mi(iShape);
            end
            
            % determine masks for y1 and y2 coordinates in whole MS

            % Compute vector masks for e.g. grad as it looks like
            % [Shape1.Dy1; Shape1.Dy2; Shape2.Dy1; Shape2.Dy2; ...
            % ShapeN.Dy1; ShapeN.Dy2]
            ptsCount = 0;
            mask_y1 = false(2*this.M,1);
            mask_y2 = false(2*this.M,1);
            
            for iShape = 1:this.nShapes
                mask_y1(ptsCount + 1 : ptsCount + Mi(iShape)) = true;
                ptsCount = ptsCount + Mi(iShape);
                mask_y2(ptsCount + 1 : ptsCount + Mi(iShape)) = true;
                ptsCount = ptsCount + Mi(iShape);
            end
            
            this.y1Mask = mask_y1;
            this.y2Mask = mask_y2;
            
        end
        
        %------------------------------------------------------------------
        % Differentiation matrices
        %------------------------------------------------------------------
        
        function ComputeDifferentiationMatrices(this)
        % Differentiation matrices are just stacked blockwise

            % Compute differentiation matrices for each shape
            for iShape = 1:this.nShapes
                this.Shapes(iShape).Shape.ComputeDifferentiationMatrix;
            end
            
            % fields to compute
            diffFields = {'Dy1','Dy2','DDy1','DDy2','Dy1Dy2','Lap','grad','div'};
            
            % loop through fields and shapes and assign in block diagonal
            % manner
            for iField = 1:length(diffFields)
                allD = [];
                for iShape = 1:this.nShapes
                    shapeD = this.Shapes(iShape).Shape.Diff.(diffFields{iField});
                    allD = blkdiag(allD, shapeD);
                end
                this.Diff.(diffFields{iField}) = allD;
            end
            
        end
        
        %------------------------------------------------------------------
        % Interpolation matrices
        %------------------------------------------------------------------

        function ComputeInterpolationMatrices(this,interpGrid)
            
            % default grid
            if(nargin<2)
                interpGrid = (-1:0.01:1)';
            end
            
            % initialisation
            InterPol = [];
            pts1 = [];
            pts2 = [];
            
            % compute interpolation Shapewise using the same interpGrid (in
            % computation space) for each Shape
            for iShape = 1:this.nShapes

                IPi = this.Shapes(iShape).Shape.ComputeInterpolationMatrix(interpGrid,interpGrid,true,true);
                InterPoli = IPi.InterPol;
                
                pts1i = IPi.pts1;
                pts2i = IPi.pts2;
                
                % convert to Cartesian points for polar shapes as pts1i,
                % pts2i are in polars for the Shape
                if strcmp(this.Shapes(iShape).Shape.polar,'polar')
                    ptsCart = GetCartPts(this.Shapes(iShape).Shape,pts1i,pts2i);
                    pts1i = ptsCart.y1_kv;
                    pts2i = ptsCart.y2_kv;
                end
                
                % MS interpolation matrix is block diagonal
                InterPol = blkdiag(InterPol,InterPoli);
                
                % MS interpolation points are a two vectors (for y1, y2)
                pts1 = [pts1;pts1i]; %#ok
                pts2 = [pts2;pts2i]; %#ok
            end
            
            % assign to MS
            this.Interp.InterPol = InterPol;
            this.Interp.pts1 = pts1;
            this.Interp.pts2 = pts2;
            
        end
        

        %------------------------------------------------------------------
        % Interpolation matrices for uniform points
        %------------------------------------------------------------------

        function InterPol = ComputeInterpolationMatricesUniform(this,opts)
            
            if(nargin<2)
                yLims = this.GetBounds;
                N1plot = 30;
                N2plot = 30;
            else
                yLims = opts.yLims;
                N1plot = opts.N1;
                N2plot = opts.N2;
            end
            
            y1Min = yLims.y1Min;
            y1Max = yLims.y1Max;
            y2Min = yLims.y2Min;
            y2Max = yLims.y2Max;
            
            dy1 = (y1Max-y1Min)/(N1plot-1);
            dy2 = (y2Max-y2Min)/(N2plot-1);

            y1_plot = (y1Min:dy1:y1Max)';
            y2_plot = (y2Min:dy2:y2Max)';
            
            pts.y1_kv = kron(y1_plot,ones(length(y2_plot),1));
            pts.y2_kv = kron(ones(length(y1_plot),1),y2_plot);
                        
            interp = ComputeInterpolationMatricesPhys(this,pts);

            this.Interp.InterPolUniform = interp.InterPol;
            this.Interp.y1Uniform = interp.pts1;
            this.Interp.y2Uniform = interp.pts2;

            InterPol = interp.InterPol;
            
        end

        %------------------------------------------------------------------
        % Interpolation matrices for physical points
        %------------------------------------------------------------------        
        
        function Interp = ComputeInterpolationMatricesPhys(this,Pts)
                        
            % points to interpolate onto
            Y1 = Pts.y1_kv;
            Y2 = Pts.y2_kv;
                        
            % initialise
            InterPol = zeros(length(Y1),length(this.Pts.y1_kv));
            
            % keep track of which point we've interpolated onto
            doneMask = false(size(Y1));
            keepMaskFull = false(size(Y1));
            
            % loop through shapes and determine which points lie in the
            % shape
            
            for iShape = 1:this.nShapes

                % get extent of shape in computational space
                InterPolShape = [];

                thisShape = this.Shapes(iShape).Shape;
                
                x1 = thisShape.Pts.x1;
                x2 = thisShape.Pts.x2;
                xTol = 10^-10;
                x1Min = min(x1) - xTol;
                x1Max = max(x1) + xTol;
                x2Min = min(x2) - xTol;
                x2Max = max(x2) + xTol;

                % transform all points from physical to computational space
                % and deal with polar coordinates
                if(strcmp(thisShape.polar,'polar'))
                    Origin = thisShape.Origin;
                    [Theta,R] = cart2pol(Y1-Origin(1),Y2-Origin(2));
                    Theta = mod(Theta,2*pi);
                    [X1,X2] = this.Shapes(iShape).Shape.CompSpace(R,Theta);
                else
                    [X1,X2] = this.Shapes(iShape).Shape.CompSpace(Y1,Y2);
                end
                
                % only keep those which came from this shape and which we
                % haven't already done
                keepMask = (X1>=x1Min) & (X1<=x1Max) & (X2>=x2Min) & (X2<=x2Max) & ~doneMask;
                keepMaskFull = keepMaskFull | keepMask;
                X1_keep = X1(keepMask);
                X2_keep = X2(keepMask);
                
                % compute interpolation matrix using shape
                nPoints = length(X1_keep);
                for iPoint = 1:nPoints
                    InterpShape = thisShape.ComputeInterpolationMatrix(X1_keep(iPoint),X2_keep(iPoint));
                    InterPolShape = [InterPolShape;InterpShape.InterPol]; %#ok
                end
                
                % combine all together
                % act on points in this shape and map to output points
                InterPol(keepMask,this.Shapes(iShape).PtsMask) = InterPolShape;
                
                % keep track of which ones points we've done (deals with
                % intersections)
                doneMask(keepMask) = true;
                
            end
            
            % throw away points not in any of the shapes
            InterPol(~keepMaskFull,:) = [];
            y1Keep = Pts.y1_kv(keepMaskFull);
            y2Keep = Pts.y2_kv(keepMaskFull);
            
            Interp.InterPol = InterPol;
            Interp.pts1 = y1Keep;
            Interp.pts2 = y2Keep;
%             Interp.pts1 = Pts.y1_kv;
%             Interp.pts2 = Pts.y2_kv;
            
        end        
        
        %------------------------------------------------------------------
        % MS Indices (complicated!)
        %------------------------------------------------------------------

        function ComputeIndices(this,normalsOverride)
            
            % compute indices for Shapes
            for iShape = 1:this.nShapes
                this.Shapes(iShape).Shape.ComputeIndices();
            end

            % list of possible edges (based on Quadrilateral and Wedge; can
            % be added to for other Shapes
            lrtbFields = {'left','right','top','bottom','radial_1','radial_2','outR','inR'};
            nFields = length(lrtbFields);
            
            %--------------------------------------------------------------
            % intersections
            %--------------------------------------------------------------

            % loop over pairs of shapes and determine which edges intersect
            
            for iShape = 1:this.nShapes
            % first shape
    
                % indices for this shape
                iInd = this.Shapes(iShape).Shape.Ind;
                
                % points for this shape
                y1i = this.Shapes(iShape).Shape.Pts.y1_kv;
                y2i = this.Shapes(iShape).Shape.Pts.y2_kv;
                
                % convert to Cartesian if necessary
                if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))
                    ptsCart = GetCartPts(this.Shapes(iShape).Shape,y1i,y2i);
                    y1i = ptsCart.y1_kv;
                    y2i = ptsCart.y2_kv;
                end
                
                for jShape = 1:this.nShapes
                % second shape
                    
                    if(iShape ~= jShape) % no self-intersections
                        
                        % mask for intersecting points between shapes i and
                        % j
                        maskij = false(this.M,1);
                        
                        % indices for this shape
                        jInd = this.Shapes(jShape).Shape.Ind;
                        
                        % points for this shape
                        y1j = this.Shapes(jShape).Shape.Pts.y1_kv;
                        y2j = this.Shapes(jShape).Shape.Pts.y2_kv;
                        
                        % convert to Cartesian if necessary
                        if(strcmp(this.Shapes(jShape).Shape.polar,'polar'))
                            ptsCart = GetCartPts(this.Shapes(jShape).Shape,y1j,y2j);
                            y1j = ptsCart.y1_kv;
                            y2j = ptsCart.y2_kv;
                        end
                        
                        % loop over pairs of fields and determine if the
                        % associated points match
                    
                        for iField = 1:nFields
                            
                            % points in the edge iField of shape iShape
                            % (may be empty)
                            if(isfield(iInd, lrtbFields{iField}))
                                mask1 = iInd.(lrtbFields{iField});
                            else
                                mask1 = [];
                            end
                            
                            y11 = y1i(mask1);
                            y21 = y2i(mask1);
                            
                            
                            for jField = 1:nFields
                                
                                % points in the edge jField of shape jShape
                                % (may be empty)

                                if(isfield(jInd,lrtbFields{jField}))
                                    mask2 = jInd.(lrtbFields{jField});
                                else
                                    mask2 = [];
                                end
                                
                                y12 = y1j(mask2);
                                y22 = y2j(mask2);
                                
                                % determine if the points match, and
                                % whether their ordering is the same or
                                % reversed (flip)
                                [matches,flip] = this.CheckIntersection(y11,y21,y12,y22);
                                
                                % store intersections between Shapes
                                if(matches)
                                    
                                    % store points in the edge of Shape
                                    % iShape (in mask1) in the correct part
                                    % of the large mask (by accessing the
                                    % part corresponding to that shape)
                                    maskij(this.Shapes(iShape).PtsMask) = mask1;
                                    this.Intersections(iShape,jShape).PtsMask = maskij;
                                    
                                    % store which side of Shape iShape is
                                    % in the intersection
                                    this.Intersections(iShape,jShape).Side = lrtbFields{iField};
                                    % and whether the sides are flipped
                                    % relative to each other
                                    this.Intersections(iShape,jShape).Flip = flip;
                                    
                                    % store corresponding corners of
                                    % intersection
                                    cornersij = false(this.M,1);
                                    cornersij(this.Shapes(iShape).PtsMask) = iInd.corners;
                                    
                                    this.Intersections(iShape,jShape).Corners = cornersij;
                                    
                                end
                                
                            end
                            
                        end
                        
                    else % iShape = jShape
                        
                        this.Intersections(iShape,jShape).PtsMask = false(this.M,1);
                        this.Intersections(iShape,jShape).Corners = false(this.M,1);
                        this.Intersections(iShape,jShape).Side = [];
                        
                    end
                    
                end
                
            end
            
            %--------------------------------------------------------------
            % end intersections
            %--------------------------------------------------------------

            %--------------------------------------------------------------
            % boundaries (doesn't store anything)
            %--------------------------------------------------------------
            % Find boundaries of MultiShape by taking the union of all of
            % the corresponding boundaries except those which are in an
            % intersection
            
            boundaries = [];
            nfields = length(lrtbFields);
            for iShape = 1:this.nShapes
                
                % Start with no edges being on the boundary of the
                % MultiShape
                boundariesTemp = {};
                
                % Add in edges if they aren't in an intersection
                for i=1:nfields
                    tempInd = this.Shapes(iShape).Shape.Ind;
                    field = lrtbFields{i};
                    % check if field exists in shape
                    if(isfield(tempInd,field))
                        % determine if the edge is already in an
                        % intersection
                        intersectionTest = false;
                        for jShape = 1:this.nShapes                            
                            if(strcmp(this.Intersections(iShape,jShape).Side,field))
                                intersectionTest = true;
                            end
                        end
                        % if not then add it to the boundary for this Shape
                        if(~intersectionTest)
                            boundariesTemp = cat(2,boundariesTemp,field);
                        end
                    end
                end
                
                % Construct contribution to the boundary from
                % Shape(iShape)
                boundaries = cat(2,boundaries,{boundariesTemp});
                
            end
            
            %--------------------------------------------------------------
            % end boundaries
            %--------------------------------------------------------------

            %--------------------------------------------------------------
            % MS edges
            %--------------------------------------------------------------
            % Loop through left, right, top and bottom and combine these
            % edges that come from each of Shape(iShape) into the
            % correseponding edge for the MultiShape

            % *** this only really makes sense for simple MS where, e.g.,
            % 'left' is well defined ***

            % loop through the edges
            for iField = 1:nFields
                
                boundary = false(this.M,1);
                
                % the current edge
                field = lrtbFields{iField};
                
                for iShape = 1:this.nShapes
                    
                    % If the edge of Shape(iShape) is in the boundary then
                    % add that to the correseponding boundary of MultiShape
                    if(any(strcmp(field,boundaries{iShape})))
                        % need to ensure that the edge mask goes into the
                        % corresponding part of the MultiShape pts by using
                        % the PtsMask for Shapes(iShape)
                        boundary(this.Shapes(iShape).PtsMask) = this.Shapes(iShape).Shape.Ind.(field);
                        
                    end
                    
                end
                
                % Save mask for MultiShape
                this.Ind.(field) = boundary;
            end
            
            % Compute entire boundary by taking intersection of all MS
            % edges
            bound = false(this.M,1);
            for iField = 1:nFields
                field = lrtbFields{iField};
                
                if(isfield(this.Ind,field))
                    bound = this.Ind.(field) | bound;
                end
                
            end
            this.Ind.bound = bound;
            
            %--------------------------------------------------------------
            % end MS edges
            %--------------------------------------------------------------

            %--------------------------------------------------------------
            % Easy indexing of subshape boundaries
            %--------------------------------------------------------------
            
            % possible edges and their normals (note corresponding order is
            % important)
            boundFields = {'left','right','top','bottom','radial_1','radial_2','outR','inR','bound'};
            normalFields = {'normalLeft','normalRight','normalTop','normalBottom','normalRadial_1',...
                'normalRadial_2','normalOutR','normalInR','normal'};
            nFields = length(boundFields);
            this.Ind.Shape = struct;
            
            % assign masks for subshape edges
            for iField = 1:nFields
                % the current edge
                boundField = boundFields{iField};
                for iShape = 1:this.nShapes
                    % need to ensure that the edge mask goes into the
                    % corresponding part of the MultiShape pts by using
                    % the PtsMask for Shapes(iShape)
                    tempInd = this.Shapes(iShape).Shape.Ind;
                    if(isfield(tempInd,boundField))
                        edge = false(this.M,1);
                        edge(this.Shapes(iShape).PtsMask) = tempInd.(boundField);
                        this.Ind.Shape(iShape).(boundField) = edge;
                    end 
                end 
            end
            
            % assign masks for subshape normals
            nFields = length(normalFields);
            for iField = 1:nFields
                
                % the current edge
                normalField = normalFields{iField};
                boundField = boundFields{iField};
                
                for iShape = 1:this.nShapes
                    % need to ensure that the normal mask goes into the
                    % corresponding part of the MultiShape pts by using
                    % the PtsMaskVec for Shapes(iShape)
                    
                    tempInd = this.Shapes(iShape).Shape.Ind;
                    if(isfield(tempInd,normalField))
                        normal = zeros(this.M,2*this.M);
                        normalCart = zeros(this.M,2*this.M);
                        
                        maski1 = this.Ind.Shape(iShape).(boundField);
                        maski2 = this.Shapes(iShape).PtsMaskVec;
                                                     
                        normal(maski1,maski2) = this.Shapes(iShape).Shape.Ind.(normalField);
                        
                        % also determine Cartesian normal for polar shapes
                        if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))    
                            normalCart(maski1,maski2) = this.Shapes(iShape).Shape.Ind.([normalField 'Cart']);
                        end
                                       
                        % delete parts of the normal that aren't on the
                        % boundary
                        normal(~maski1,:) = [];
                        this.Ind.Shape(iShape).(normalField) = sparse(normal);
                        
                        % same for Cartesian normal
                        if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))    
                            normalCart(~maski1,:) = [];
                            this.Ind.Shape(iShape).([normalField 'Cart']) = sparse(normalCart);
                        end

                    end
                    
                    
                end
                
            end

            %--------------------------------------------------------------
            % end Shapewise indexing
            %--------------------------------------------------------------

            %--------------------------------------------------------------
            % masks for intersections
            %--------------------------------------------------------------

            % Construct intersection mask by looping through Shapes
            intersections = false(this.M,1);
            for iShape = 1:this.nShapes
                for jShape = 1:this.nShapes
                    intersections(this.Intersections(iShape,jShape).PtsMask) = true;
                end
            end
            this.Ind.intersections = intersections;

            %--------------------------------------------------------------
            % end intersection masks
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % normals for MS
            %--------------------------------------------------------------
            % Construct normals for the MultiShape

            % loop over all fields and shapes
            nFields = length(boundFields);
            
            for iField = 1:nFields
                
                allNormals = [];
                allNormalsCart = [];
                
                for iShape = 1:this.nShapes
                    
                    % number of points in the Shape
                    N1i = this.Shapes(iShape).Shape.Pts.N1;
                    N2i = this.Shapes(iShape).Shape.Pts.N2;
                    Mi = N1i*N2i;
                    
                    tempInd = this.Shapes(iShape).Shape.Ind;
                    
                    % check if the field exists in the Shape
                    if(isfield(tempInd,boundFields{iField}))

                        % get boundary from Shape
                        bound = this.Shapes(iShape).Shape.Ind.(boundFields{iField});
                        % and normal from MS (which is copied from Shape)
                        normal = this.Ind.Shape(iShape).(normalFields{iField});
                        % and the Cartesian normal for polar shapes
                        if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))   
                            normalCart = this.Ind.Shape(iShape).([normalFields{iField} 'Cart']);
                        end
                        
                        % assign normals
                        temp = zeros(Mi,2*Mi);
                        temp(bound,:) = normal(:,this.Shapes(iShape).PtsMaskVec);
                        
                        % and the Cartesian version
                        if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))    
                            tempCart = zeros(Mi,2*Mi);
                            tempCart(bound,:) = normalCart(:,this.Shapes(iShape).PtsMaskVec);
                        end
                    else  % default to nothing
                        temp = zeros(Mi,2*Mi);
                        tempCart = zeros(Mi,2*Mi);
                    end

                    % stack normals together in block diagonal manner
                    allNormals = blkdiag(allNormals,temp);
                    % do the same for the Cartesian version
                    if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))    
                        allNormalsCart = blkdiag(allNormalsCart,tempCart);
                    else % for non-polar case can just use the normal version
                        allNormalsCart = blkdiag(allNormalsCart,temp);
                    end
                end

                % delete anything not in the boundary
                allNormals(~this.Ind.(boundFields{iField}),:) = []; %#ok
                this.Ind.(normalFields{iField}) = sparse(allNormals);
                
                % and for Cartesian case
                allNormalsCart(~this.Ind.(boundFields{iField}),:) = []; %#ok
                this.Ind.([normalFields{iField} 'Cart']) = sparse(allNormalsCart);

            end
            
            % Fix normals at the intersections of shapes and the boundary
            
            % use Cartesian normals
            normal = this.Ind.normalCart;
            
            % find repeated points on the boundary (corners of intersecting
            % shapes)
            tol = 10^-10;
            boundPts = [this.Pts.y1_kv(this.Ind.bound),this.Pts.y2_kv(this.Ind.bound)];
            %[~,uniqueIdx] = unique(boundPts,'rows');
            % find unique points in the boundary up to a tolerance
            [uniqueBoundPts,uniqueIdx,IC] = uniquetol(boundPts,tol,'ByRows',true);
            % reconstruct bound points up to tolerance (i.e. replace those
            % within the tolerance by the 'original')
            boundPts = uniqueBoundPts(IC,:);
            % determine duplicates
            duplicateIdx = setdiff(1:size(boundPts, 1), uniqueIdx);
            duplicates = boundPts(duplicateIdx,:);
            
            if(~isempty(normalsOverride))
                nOverride = size(normalsOverride.pts,1);
                overridePts = normalsOverride.pts;
                overrideNormals = normalsOverride.normals;
            else
                overridePts = [];
                overrideNormals = [];
                nOverride = 0;
            end
            
            % loop through the duplicates and average the normal
            originalPts = [];
            originalNormals = [];
            for iDup = 1:size(duplicates,1)
                % allow normals to be overridden by input
                normalTemp = [];
                
                for iOverride = 1:nOverride                    
                    if(this.myisequal(overridePts(iOverride,:),duplicates(iDup,:)))
                        normalTemp = overrideNormals(iOverride,:);                       
                    end
                   
                end

                % determine number of matching points
                temp = ( boundPts == duplicates(iDup,:) );
                matches = find(temp(:,1) & temp(:,2));
                nMatches = length(matches);
                
                if(isempty(normalTemp))
                
                    normalTemp = [0,0];

                    % extract normal vectors and sum
                    for iMatch = 1:nMatches
                        normalMatch = normal(matches(iMatch),:);
                        normalMatch(normalMatch==0) = [];
                        normalTemp = normalTemp + normalMatch;
                    end

                    % normalise normal vector
                    normalTempNorm = sqrt( normalTemp(:,1).^2 + normalTemp(:,2).^2 );
                    normalTempNorm2 = [normalTempNorm normalTempNorm];
                    normalTemp = normalTemp./normalTempNorm2;

                end
                 % reassign the same vector to all the relevant points
                for iMatch = 1:nMatches
                    normalMatch = normal(matches(iMatch),:);
                    % saving duplicate normals and corresponding points
                    originalNormals = [originalNormals; normalMatch];
                    originalPts = [originalPts; duplicates(iDup,:)];  
                    
                    normalMatch(normalMatch~=0) = normalTemp;
                    normal(matches(iMatch),:) = normalMatch;
                    normal(matches(iMatch),:);
                end

            end

            % put back in to normalCart
            this.Ind.normalCart = normal;
            
            % put shapewise normal (polar/cart) into Ind.normal
            for iPt = 1:size(normal,1)
                tempNormal = normal(iPt,:)';
                tempNormal = this.ConvertVectorToShapes(tempNormal);
                normal(iPt,:) = tempNormal';
            end
            this.Ind.normal = normal;
            
            this.Ind.originalPts = originalPts;
            this.Ind.originalNormals = originalNormals;
            this.Ind.overridePts = overridePts;
            this.Ind.overrideNormals = overrideNormals;
            
            %--------------------------------------------------------------
            % end MS normals
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % normals for MS intersections
            %--------------------------------------------------------------
           
            % loop through all pairs of shapes
            for iShape = 1:this.nShapes
                for jShape = 1:this.nShapes
                    
                    if(iShape ~= jShape) % ignore self intersections
                        
                        % determine if there is an intersection
                        sidei = this.Intersections(iShape,jShape).Side;
                        % and where that side is in the list
                        iPos = find(strcmp(sidei,boundFields));
                        
                        if(~isempty(iPos))
                            
                            normal = zeros(this.M,2*this.M);
                            normalCart = zeros(this.M,2*this.M);
                            
                            % mask for scalars
                            maski1 = this.Intersections(iShape,jShape).PtsMask; 
                            % mask for vectors
                            maski2 = this.Shapes(iShape).PtsMaskVec;
    
                            % normal corresponding to intersecting edge
                            normali = this.Shapes(iShape).Shape.Ind.(normalFields{iPos}); 
                            
                            % assign to normal and remove points not in
                            % edge
                            normal(maski1,maski2) = normali;
                            normal(~maski1,:) = [];

                            % assign to MS
                            this.Intersections(iShape,jShape).Normal = sparse(normal);
                            
                            % do the same for Cartesian normal for polar
                            % case
                            if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))
                                normaliCart = this.Shapes(iShape).Shape.Ind.([normalFields{iPos} 'Cart']); 
                                normalCart(maski1,maski2) = normaliCart;
                                normalCart(~maski1,:) = [];
                                this.Intersections(iShape,jShape).NormalCart = sparse(normalCart);
                            end
             
                        end
                        
                    end
                    
                end
            end
            
        end
        
        %------------------------------------------------------------------
        % Integration vectors
        %------------------------------------------------------------------
        
        function ComputeIntegrationVectors(this)
            
            this.Int = [];
            
            % compute Shape integration vector and stack into row vector
            for iShape = 1:this.nShapes
                Inti = this.Shapes(iShape).Shape.ComputeIntegrationVector;
                this.Int = [this.Int, Inti];
            end
        end

        %------------------------------------------------------------------
        % Convolution matrix
        %------------------------------------------------------------------
        % this is essentially copied from Shape

        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end
            
            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;
            
            if(useDistance)
                fPTemp = f(GetDistance(this,y1,y2));
            else
                fPTemp = f(y1,y2);
            end
            
            fDim = size(fPTemp);
            
            nElts = prod(fDim(2:end));
            
            IntT = (this.Int).';  % M x 1
            
            IntT = IntT(:,ones(1,nElts)); % M x nElts
            IntT = reshape(IntT,fDim);    % size(f)
            
            M_conv = zeros([this.M,this.M,fDim(2:end)]);
            
            Mmask = repmat({':'},[1,fDim]);
            
            for i=1:this.M
                if(useDistance)
                    fP          = f(GetDistance(this,y1(i) - y1,y2(i) - y2));
                else
                    fP          = f(y1(i) - y1,y2(i) - y2);
                end
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            M_conv(isnan(M_conv)) = 0;
            
            if((nargin >= 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
            
        end

        %------------------------------------------------------------------
        % Euclidean distance
        %------------------------------------------------------------------
        % this is essentially copied from Shape

        function d = GetDistance(this,y1,y2) %#ok
            d = sqrt( y1.^2 + y2.^2 );
        end

        %------------------------------------------------------------------
        % Apply intersection boundary conditions
        %------------------------------------------------------------------

        function dydt = ApplyIntersectionBCs(this,dydt,data,BCs)
        % note this modifies the given dydt used in a DAE solver

            % loop through pairs of shapes
            for iShape = 1:this.nShapes
                for jShape = iShape+1:this.nShapes

                    % if they intersect apply the given BC function
                    if(~isempty(BCs(iShape,jShape).function))                        
                        BCfn = str2func(BCs(iShape,jShape).function);
                        [BC,mask] = BCfn(data,this.Intersections,iShape,jShape);
                        dydt(mask) = BC(mask);
                    end
                end
            end
            
        end

        %------------------------------------------------------------------
        % Make and split MS vectors given the two coordinates
        %------------------------------------------------------------------
        
        function v = MakeVector(this,v1,v2)
            v = zeros(2*length(v1),1);
            v(this.y1Mask) = v1;
            v(this.y2Mask) = v2;
        end
        
        function [v1,v2] = SplitVector(this,v)
            v1 = v(this.y1Mask);
            v2 = v(this.y2Mask);
        end

        %------------------------------------------------------------------
        % Convert a vector to/from each Shape's coordinate system
        %------------------------------------------------------------------        

        function vCartPolar = ConvertVectorToShapes(this,vCart)

            % get vector in coordinates
            [vCart1,vCart2] = SplitVector(this,vCart);
            
            % copy into output vector
            vCartPolar1 = vCart1;
            vCartPolar2 = vCart2;
            
            % convert the parts on polar shapes to polar coordinates
            for iShape = 1:this.nShapes
                
                if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))
                   
                    % theta for shape
                    mask = this.Shapes(iShape).PtsMask;
                    theta = this.Shapes(iShape).Shape.Pts.y2_kv;
                    
                    % Cartesian vector components
                    v1 = vCart1(mask);
                    v2 = vCart2(mask);

                    % polar vector components
                    vr = cos(theta).*v1 + sin(theta).*v2;
                    vtheta = - sin(theta).*v1 + cos(theta).*v2;
                    vCartPolar1(mask) = vr;
                    vCartPolar2(mask) = vtheta;
                    
                end
            end
            
            % reassemble into vector
            vCartPolar = MakeVector(this,vCartPolar1,vCartPolar2);
            
        end

        function vCart = ConvertVectorFromShapes(this,vCartPolar)

            % get vector in coordinates
            [vCartPolar1,vCartPolar2] = SplitVector(this,vCartPolar);
            
            % copy into output vector
            vCart1 = vCartPolar1;
            vCart2 = vCartPolar2;
            
            % convert the parts on polar shapes to polar coordinates
            for iShape = 1:this.nShapes
                
                if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))
                   
                    % theta for shape
                    mask = this.Shapes(iShape).PtsMask;
                    theta = this.Shapes(iShape).Shape.Pts.y2_kv;
                    
                    % Polar vector components
                    vr = vCartPolar1(mask);
                    vtheta = vCartPolar2(mask);

                    % Cartesian vector components
                    v1 = cos(theta).*vr - sin(theta).*vtheta;
                    v2 = sin(theta).*vr + cos(theta).*vtheta;
                    vCart1(mask) = v1;
                    vCart2(mask) = v2;
                    
                end
            end
            
            % reassemble into vector
            vCart = MakeVector(this,vCart1,vCart2);
            
        end
        
        
        
        %------------------------------------------------------------------
        % Plotting
        %------------------------------------------------------------------

        
        %------------------------------------------------------------------
        % Plot grid
        %------------------------------------------------------------------

        function PlotGrid(this)
            % plot grid for each Shape
            for iShape = 1:this.nShapes
                thisShape = this.Shapes(iShape).Shape;
                thisShape.PlotGrid;
                hold on
            end
        end

        
        %------------------------------------------------------------------
        % Plot MS bound
        %------------------------------------------------------------------

        function PlotBound(this,boundList,displayNames)
        % boundList is a list of edges to plot, default is to just plot
        % bound

            if(nargin==1)
                boundList = {'bound'};
                displayNames = {'Boundary'};
            end
            
            if(nargin==2)
                displayNames = boundList;
            end
            
            colourList = {'r','b','m','c'};
            
            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;
            
            nFields = length(boundList);
            
            % loop through fields to plot and plot markers on the edges
            for iField = 1:nFields
                mask = this.Ind.(boundList{iField});
                p1 = scatter(y1(mask),y2(mask), ...
                    'o','MarkerEdgeColor',colourList{iField}, ...
                    'MarkerFaceColor',colourList{iField}, ...
                    'DisplayName',displayNames{iField});
                hold on
            end
        end
        
        %------------------------------------------------------------------
        % Plot MS intersections
        %------------------------------------------------------------------

        function PlotIntersections(this,plotIJ,opts)
        % plotIJ is logical(nShapes,nShapes) for which to plot; default is
        % all of them

            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;
            
            if(nargin<2 || isempty(plotIJ))
                plotIJ = true(this.nShapes);
            end
            
            if(nargin<3)
                opts = [];
            end
            
            % get options
            if(isfield(opts,'colour'))
                colour = opts.colour;
            else
                colour = 'b';
            end
            
            if(isfield(opts,'linewidth'))
                linewidth = opts.linewidth;
            else
                linewidth = 2;
            end
            
            if(isfield(opts,'plotPoints'))
                plotPoints = opts.plotPoints;
            else
                plotPoints = false;
            end
            
            firstPlot = true;
            % loop over all possible intersections
            for iShape = 1:this.nShapes
                for jShape = iShape+1:this.nShapes
                    if(plotIJ(iShape,jShape))
                        PtsMask = this.Intersections(iShape,jShape).PtsMask;
                        % scatter plot
                        if(plotPoints)
                            scatter(y1(PtsMask),y2(PtsMask),'x',colour);
                            hold on
                        end
                        % line plot
                        if(firstPlot)
                            plot(y1(PtsMask),y2(PtsMask),'Color',colour,'Linewidth',linewidth, ...
                                 'DisplayName','Intersections');
                            firstPlot = false;
                        else
                            plot(y1(PtsMask),y2(PtsMask),'Color',colour,'Linewidth',linewidth, ...
                                 'HandleVisibility','off');
                        end
                    end
                end
            end
            
            % fix axes
            this.SetAxes;
            
        end
        
        %------------------------------------------------------------------
        % Plot Shape normals
        %------------------------------------------------------------------

        function PlotNormalsShapes(this)
            
            % loop through shapes
            for iShape = 1:this.nShapes

                % shape points
                y1 = this.Shapes(iShape).Shape.Pts.y1_kv;
                y2 = this.Shapes(iShape).Shape.Pts.y2_kv;
                
                % Cartesian points for polar shape
                if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))
                    ptsCart = this.Shapes(iShape).Shape.GetCartPts(y1,y2);
                    y1 = ptsCart.y1_kv;
                    y2 = ptsCart.y2_kv;
                end
                
                % Cartesian directions
                dir1 = [ones(size(y1)); zeros(size(y2))];
                dir2 = [zeros(size(y1)); ones(size(y2))];
                
                % normal
                if(strcmp(this.Shapes(iShape).Shape.polar,'polar'))
                    normal = this.Shapes(iShape).Shape.Ind.normalCart;
                else
                    normal = this.Shapes(iShape).Shape.Ind.normal;
                end
            
                % boundary
                bound = this.Shapes(iShape).Shape.Ind.bound;
                
                % construct and normalise vector
                v1 = normal * dir1;
                v2 = normal * dir2;
                
                vnorm = sqrt(v1.^2 + v2.^2);
                v1 = v1./vnorm;
                v2 = v2./vnorm;
                
                % cut down to points on the corresponding boundary
                y1 = y1(bound);
                y2 = y2(bound);
                
                % plot arrows with forced normalisation
                quiver(y1, y2, v1, v2, 'AutoScale','off');
                
                hold on
                
            end
            
            % set axes so the arrows are scaled correctly
            axis equal
            
        end
        

        %------------------------------------------------------------------
        % Plot MS normals
        %------------------------------------------------------------------

        function PlotNormals(this,boundList)

            % list of possible boundaries and their normals
            boundListFull = {'left','right','top','bottom','radial_1','radial_2','outR','inR','bound'};
            normalListFull = {'normalLeft','normalRight','normalTop', ...
                'normalBottom','normalRadial_1',...
                'normalRadial_2','normalOutR','normalInR','normal'};
            
            % defaults
            if(nargin==1)
                boundList = {'bound'};
            end
            
            colourList = {'r','b','m','c','k'};
            

            % points
            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;
            
            % Cartesian unit vectors
            dir1 = zeros(2*this.M,1);
            dir2 = zeros(2*this.M,1);
            dir1(this.y1Mask) = 1;
            dir2(this.y2Mask) = 1;
            
            % loop over given fields and extract points and normals for MS
            nFields = length(boundList);
            
            for iField = 1:nFields 
                iPos = find(strcmp(boundList{iField},boundListFull));
                if(~isempty(iPos))
                    % select Cartesian version of normals
                    if(isfield(this.Ind,[normalListFull{iPos} 'Cart']))
                        normal = this.Ind.([normalListFull{iPos} 'Cart']);
                    else
                        normal = this.Ind.(normalListFull{iPos});
                    end
                    % points
                    bound = this.Ind.(boundListFull{iPos});
                end
                
                % construct normalised vectors
                v1 = normal * dir1;
                v2 = normal * dir2;
                vnorm = sqrt(v1.^2 + v2.^2);
                v1 = v1./vnorm;
                v2 = v2./vnorm;

                % and plot them on the boundary
                y1Plot = y1(bound);
                y2Plot = y2(bound);

                quiver(y1Plot, y2Plot, v1, v2, 'AutoScale','off', ...
                       'Color',colourList{mod(iPos,5)+1},'LineWidth',1.5, ...
                       'DisplayName','Normals');
                hold on
             end
            axis equal
        end
        
        %------------------------------------------------------------------
        % Plot MS normals to intersections
        %------------------------------------------------------------------

        function PlotNormalsIntersections(this)
            
            % points
            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;
            
            % unit vectors
            dir1 = zeros(2*this.M,1);
            dir2 = zeros(2*this.M,1);
            dir1(this.y1Mask) = 1;
            dir2(this.y2Mask) = 1;
            
            % loop over all possible pairs of shapes
            for iShape = 1:this.nShapes
                for jShape = 1:this.nShapes
                    
                    if(any(this.Intersections(iShape,jShape).PtsMask))  
                        % find Cartesian normal for interesctions
                        if(isfield(this.Intersections(iShape,jShape),'NormalCart'))
                            normal = this.Intersections(iShape,jShape).NormalCart;
                        else
                            normal = this.Intersections(iShape,jShape).Normal;
                        end
                            
                        % and corresponding points
                        bound = this.Intersections(iShape,jShape).PtsMask;
                            
                        % normalise vectors
                        v1 = normal * dir1;
                        v2 = normal * dir2;
                        
                        vnorm = sqrt(v1.^2 + v2.^2);
                        v1 = v1./vnorm;
                        v2 = v2./vnorm;
                        
                        % plot vectors at bound points
                        y1Plot = y1(bound);
                        y2Plot = y2(bound);
                        
                        quiver(y1Plot, y2Plot, v1, v2, 'AutoScale','off','LineWidth',1.5);
                        
                        hold on
                    end
                    
                end
                
            end
            
            axis equal
            
        end
        
        %------------------------------------------------------------------
        % Plot function on MS
        %------------------------------------------------------------------

        function Plot(this,V,opts,optsDetail)
            
            if(nargin<4)
                optsDetail = [];
            end
            if(nargin<3)
                opts = {};
            end
            
            optsDetail.sigma = false;
            optsDetail.edgecolor = 'none';
            
            % plot Shapewise
            for iShape = 1:this.nShapes
                % cut down to points in the Shape
                Vi = V(this.Shapes(iShape).PtsMask);
                % plot the result
                this.Shapes(iShape).Shape.plot(Vi,opts,optsDetail);
                hold on
            end
            
            this.SetAxes;
            
        end

        %------------------------------------------------------------------
        % Plot vector field on MS
        %------------------------------------------------------------------
        
        function PlotVectorField(this,v,opts)
            
            if(nargin<3)
                opts = [];
            end
            
            if(isfield(opts,'colour'))
                colour = opts.colour;
            else
                colour = 'b';
            end
            
            if(isfield(opts,'lineWidth'))
                lineWidth = opts.lineWidth;
            else
                lineWidth = 2;
            end
            
            % construct a uniform grid to plot at
            if(~isfield(this.Interp,'InterPolUniform'))
                if(~isempty(opts))
                    IPU = this.ComputeInterpolationMatricesUniform(opts);
                else
                    IPU = this.ComputeInterpolationMatricesUniform();
                end
            else
                IPU = this.Interp.InterPolUniform;
            end
            
            %PtsU = this.Pts.Uniform;
            y1U = this.Interp.y1Uniform;
            y2U = this.Interp.y2Uniform;
            
            % determine components of vector field
            v = this.ConvertVectorFromShapes(v); % need everything in Cartesians
            [v1,v2] = this.SplitVector(v);
            
            % interpolate onto plotting points
            v1U = IPU*v1;
            v2U = IPU*v2;

            % plot
            quiver(y1U,y2U,v1U,v2U,'color',colour,'linewidth',lineWidth);

            this.SetAxes;
            
        end
        
        %------------------------------------------------------------------
        % Default axes limits
        %------------------------------------------------------------------

        function SetAxes(this,ha)
            
            if (nargin < 2)
                ha = gca;
            end
            
            % get bounds
            xMin = this.Pts.y1Min;
            xMax = this.Pts.y1Max;
            yMin = this.Pts.y2Min;
            yMax = this.Pts.y2Max;
            
            % restrict to bounds
            xlim([xMin, xMax]);
            ylim([yMin, yMax]);
            
            % set aspect ratio
            pbaspect(ha,[(xMax-xMin) (yMax-yMin) 1/2*min((xMax-xMin),(yMax-yMin))]);
        end
        
        %------------------------------------------------------------------
        % Construct uniform grid over all of MS
        %------------------------------------------------------------------

        function yLims = GetBounds(this)
           
            % get bounding rectangle
            yLims.y1Min = this.Pts.y1Min;
            yLims.y1Max = this.Pts.y1Max;
            yLims.y2Min = this.Pts.y2Min;
            yLims.y2Max = this.Pts.y2Max;

        end

        %------------------------------------------------------------------
        % Check if two vectors are equal with lower tolerance than isequal
        %------------------------------------------------------------------
        
        function matches = myisequal(this,y1,y2)
            
            matches = 0;
            
            if(size(y1) == size(y2))
                d = abs(y1-y2);
                if max(d)< 10^(-14)
                    matches = 1;
                end
            end
        end

        %------------------------------------------------------------------
        % Determine if two sets of points are the same, possibly in reverse
        % order
        %------------------------------------------------------------------
        
        function [matches,flip] = CheckIntersection(this,y11,y12,y21,y22)
            
            matches = false;
            flip = false;
            
            % can't match if empty
            if(isempty(y11) || isempty(y12) || isempty(y21) || isempty(y22))
                return
            end
            
            % construct pairs of points
            y1 = [y11 y12];
            y2 = [y21 y22];
            
            % check if they match
            matches1 = this.myisequal(y1,y2);
            if(matches1)
                matches = true;
            else
                % check if they match with the order of one reversed
                matches2 = this.myisequal(y1,flipud(y2)); 
                if(matches2)
                    matches = true;
                    flip = true;
                end
            end
        end
        
    end % end methods
    
end