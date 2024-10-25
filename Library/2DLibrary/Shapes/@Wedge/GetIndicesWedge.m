function I  = GetIndicesWedge(this)
%************************************************************************
%I  = GetIndicesBox(x1,x2)
%INPUT:
% x1 - grid points in computational space in 1st variable, x1 in [-1 1]
% x2 - grid points in computational space in 2nd variable, x2 in [-1 1]
%OUTPUT: structure with ...
% - indices of 'radial1', 'radial2', 'inner', 'outer' boundaries
% - bound = includes radial_1, radial_2, outer and inner boundaries
% - normal vectors: normalRadial_1, normalRadial_2, normalouter, normalinner
% - normal = normal vectors for full boundary
% - corners
% - indices of 'finite', 'infinite' boundaries
% - normal vectors: normalFinite, normalInfinite
%************************************************************************
%         N1  = this.Pts.N1;
%         N2  = this.Pts.N2;

x1_kv = this.Pts.x1_kv; % r in comp space
x2_kv = this.Pts.x2_kv; % theta in comp space
y1_kv = this.Pts.y1_kv; % r in phys space
y2_kv = this.Pts.y2_kv; % theta in phys space

% make logicals first
outer   = (x1_kv == max(x1_kv)); % out radius = const
inner   = (x1_kv == min(x1_kv)); % in radius = const
radial1 = (x2_kv == min(x2_kv)); % theta1 = const
radial2 = (x2_kv == max(x2_kv)); % theta2 = const

bound   = (radial1 | outer | radial2 | inner);

radial1Finite  = any(isfinite(y2_kv(radial1)));
radial2Finite  = any(isfinite(y2_kv(radial2)));
outerFinite    = any(isfinite(y1_kv(outer)));
innerFinite    = any(isfinite(y2_kv(inner)));

finite1 = ( (radial1Finite & radial1) | (radial2Finite & radial2) );
finite2 = ( (outerFinite & outer) | (innerFinite & inner) );
finite  = (finite1 | finite2);

infinite = ( (~radial1Finite & radial1) | (~radial2Finite & radial2) ...
    | (~outerFinite & outer) | (~innerFinite & inner) );

corners = ((radial1 & outer) | (outer & radial2) | (radial2 & inner) | (inner & radial1));

%------------------------------------------------------------------------
% construct normals for plotting
%------------------------------------------------------------------------

Z = zeros(length(x1_kv));

% normals in r (y1) direction
nx1radial2 = Z;
nx1radial1 = Z;
nx1outer   = Z;
nx1inner   = Z;

% normals in theta (y2) direction
nx2radial2 = Z;
nx2radial1 = Z;
nx2outer   = Z;
nx2inner   = Z;

outerradial2 = (outer & radial2);
outerradial1 = (outer & radial1);
innerradial2 = (inner & radial2);
innerradial1 = (inner & radial1);


% Outer normal

n2Temp = y2_kv(outer); % thetas
n1Temp = ones(size(n2Temp)); % rs
%n1Temp = ones(this.N2,1); % rs

nx1outer(outer,outer)     = diag(n1Temp);
nx2outer(outer,outer)     = diag(n2Temp);

% convert to cartesian
normalOutRPol  = sparse( [nx1outer(outer,:)  nx2outer(outer,:)] );
normalOutRCart = this.GetCartVec(normalOutRPol);

% Radial2 normal

n2Temp = y2_kv(radial2); % thetas
n1Temp = ones(size(n2Temp)); %rs
%n1Temp = ones(this.N1,1); % rs

nx1radial2(radial2,radial2) = diag(n1Temp);

interiorradial2 = true(this.N1,1);

nx2radial2(radial2,radial2)  = diag(n2Temp)...
    +diag(pi/2*(interiorradial2));

% sparse(nx2radial2(radial2,radial2))

% if(this.Invert)
%     nx2radial2(radial2,radial2) = nx2radial2(radial2,radial2)-diag(pi/2*interiorradial2);
% end
% sparse(nx2radial2(radial2,radial2))

% convert to cartesian
normalRadial_2Pol = sparse( [nx1radial2(radial2,:) nx2radial2(radial2,:)] );
normalRadial_2Cart = this.GetCartVec(normalRadial_2Pol);

% if(this.Invert)
% normalRadial_2Cart = normalRadial_2Cart;
% end


% Inner normal

n2Temp = y2_kv(inner); % thetas
n1Temp = ones(size(n2Temp)); % rs

nx1inner(inner,inner)     = diag(n1Temp);

interiorinner = true(this.N2,1);

nx2inner(inner,inner)     = diag(n2Temp)...
    +diag(pi*(interiorinner));

% convert to cartesian
normalInRPol = sparse( [nx1inner(inner,:)  nx2inner(inner,:)] );
normalInRCart = this.GetCartVec(normalInRPol);

% Radial1 normal

n2Temp = y2_kv(radial1); % thetas
n1Temp = ones(size(n2Temp)); % rs

nx1radial1(radial1,radial1)     = diag(n1Temp);

interiorradial1 = true(this.N1,1);

nx2radial1(radial1,radial1)     = diag(n2Temp)...
    +diag(3*pi/2*(interiorradial1));


% if(this.Invert)
% nx2radial1(radial1,radial1) = nx2radial1(radial1,radial1)+diag(pi/2*(interiorradial1));
% nx2radial1(radial1,radial1)
% end

normalRadial_1Pol = sparse( [nx1radial1(radial1,:) nx2radial1(radial1,:)] );
normalRadial_1Cart = this.GetCartVec(normalRadial_1Pol);


%%% corners

%         innerradial1Ang = (nx2inner(innerradial1,innerradial1) + nx2radial1(innerradial1,innerradial1))/2
%         outerradial1Ang = (nx2radial1(outerradial1,outerradial1) + nx2outer(outerradial1,outerradial1))/2
%         outerradial2Ang = (nx2radial2(outerradial2,outerradial2) + nx2outer(outerradial2,outerradial2))/2;
%         innerradial2Ang = (nx2radial2(innerradial2,innerradial2) + nx2inner(innerradial2,innerradial2))/2
%
%         nx2inner(innerradial1,innerradial1) = innerradial1Ang;
%         nx2radial1(innerradial1,innerradial1) = innerradial1Ang;
%         nx2radial1(outerradial1,outerradial1) = outerradial1Ang;
%         nx2outer(outerradial1,outerradial1) = outerradial1Ang;
%         nx2radial2(outerradial2,outerradial2) = outerradial2Ang;
%         nx2outer(outerradial2,outerradial2) = outerradial2Ang;
%         nx2radial2(innerradial2,innerradial2) = innerradial2Ang;
%         nx2inner(innerradial2,innerradial2) = innerradial2Ang;


psiOuterRadial1 = cornerAngle(nx2radial1(outerradial1,outerradial1), nx2outer(outerradial1,outerradial1));
psiOuterRadial2 = cornerAngle(nx2radial2(outerradial2,outerradial2), nx2outer(outerradial2,outerradial2));
psiInnerRadial2 = cornerAngle(nx2radial2(innerradial2,innerradial2), nx2inner(innerradial2,innerradial2));
psiInnerRadial1 = cornerAngle(nx2radial1(innerradial1,innerradial1), nx2inner(innerradial1,innerradial1));

nThetas = nx2outer(bound,:) + nx2inner(bound,:) + nx2radial2(bound,:) + nx2radial1(bound,:);
%         nThetas = nx2outer(bound,:) | nx2inner(bound,:) | nx2radial2(bound,:) | nx2radial1(bound,:);

[~, outerradial1Corner] = ismember(find(outerradial1),find(bound));
[~, outerradial2Corner] = ismember(find(outerradial2),find(bound));
[~, innerradial2Corner] = ismember(find(innerradial2),find(bound));
[~, innerradial1Corner] = ismember(find(innerradial1),find(bound));

nThetas(outerradial1Corner,:) = psiOuterRadial1;
nThetas(outerradial2Corner,:) = psiOuterRadial2;
nThetas(innerradial2Corner,:) = psiInnerRadial2;
nThetas(innerradial1Corner,:) = psiInnerRadial1;

nRadii = nx1outer(bound,:) | nx1inner(bound,:) | nx1radial2(bound,:) | nx1radial1(bound,:);

normalPol = sparse([nRadii, nThetas]);

normalCart = this.GetCartVec(normalPol);


%------------------------------------------------------------------------
% construct normals for polar computations (these are the useful ones)
%------------------------------------------------------------------------
Z = zeros(length(x1_kv));

nxThetaRadial2 = Z;
nxThetaRadial1 = Z;

nxrOuter = Z;
nxrInner = Z;

nxThetaRadial2(radial2,radial2) = speye(sum(radial2));
nxThetaRadial1(radial1,radial1) = -speye(sum(radial1));

nxrOuter(outer,outer) = speye(sum(outer));
nxrInner(inner,inner) = -speye(sum(inner));
%------------------------------------------------------------------------
% construct normals finite and infinite
%------------------------------------------------------------------------

nf1 = radial2Finite*nx1radial2(finite,:) + radial1Finite*nx1radial1(finite,:) + ...
    outerFinite*nx1outer(finite,:)   + innerFinite*nx1inner(finite,:);
nf2 = radial2Finite*nx2radial2(finite,:) + radial1Finite*nx2radial1(finite,:) + ...
    outerFinite*nx2outer(finite,:)   + innerFinite*nx2inner(finite,:);

normalFinite = sparse([nf1 nf2]);

normalFinite1 = sparse( nf1 );
normalFinite2 = sparse( nf2 );

normalInfinite = sparse([Z Z]);

% normal = sparse([nRadii, nThetas]);

%------------------------------------------------------------------------
% Invert if required
%------------------------------------------------------------------------

% if(this.Invert)
%     
% %     radial2temp = radial2;
% %     radial1temp = radial1;
% %     innerTemp = inner;
% %     outerTemp = outer;
% %     
%     normalRadial_2CartTemp = normalRadial_2Cart;
%     normalRadial_1CartTemp = normalRadial_1Cart;
%     normalOutRCartTemp = normalOutRCart;
%     normalInRCartTemp = normalInRCart;
%     
% %     radial1 = radial2temp;
% %     radial2 = radial1temp;
% %     inner = outerTemp;
% %     outer = innerTemp;
% %     
%     normalRadial_1Cart = normalRadial_2CartTemp;
%     normalRadial_2Cart = normalRadial_1CartTemp;
%     normalInRCart = normalOutRCartTemp;
%     normalOutRCart = normalInRCartTemp;
% end

normalCart = normalCart - [this.Origin(1)*true(size(normalCart(:,1:end/2))), this.Origin(2)*true(size(normalCart(:,end/2+1:end)))];

%------------------------------------------------------------------------
% Put everything into a structure
%------------------------------------------------------------------------

I = struct('radial_2',radial2,'radial_1',radial1,'inR',inner,'outR',outer,...
        'bound',bound,...
        'normalRadial_2',sparse([Z(radial2,:)  nxThetaRadial2(radial2,:)]),...
        'normalRadial_1',sparse([Z(radial1,:)  nxThetaRadial1(radial1,:)]),...
        'normalOutR',sparse([nxrOuter(outer,:)  Z(outer,:)]),...
        'normalInR',sparse([nxrInner(inner,:)  Z(inner,:)]),...
        'normal',sparse([(nxrOuter(bound,:)+nxrInner(bound,:)) (nxThetaRadial1(bound,:) + nxThetaRadial2(bound,:))]),...
        'normalRadial_2Cart', normalRadial_2Cart,...
        'normalRadial_1Cart', normalRadial_1Cart,...
        'normalOutRCart', normalOutRCart,...
        'normalInRCart', normalInRCart,...
        'normalCart', normalCart,...
        'corners',corners, ...
        'finite',finite,'infinite',infinite, ...
        'normalFinite',normalFinite,'normalInfinite',normalInfinite,...
        'finite1',finite1,'finite2',finite2,...
        'normalFinite1',normalFinite1,'normalFinite2',normalFinite2);


    function psi = cornerAngle(theta1,theta2)
        
        phi = theta2 - theta1;
        if(abs(phi)<=pi)
            psi = theta1 + phi/2; % simple case
        elseif pi<phi && phi<3/2*pi
            psi = theta1 + (2*pi-phi)/2; % flipped case 3rd quadrant
        else
            psi = theta1-(2*pi-phi)/2;  % flipped case 4th quadrant
        end
        psi = mod(psi,2*pi); % probably not needed
    end



end