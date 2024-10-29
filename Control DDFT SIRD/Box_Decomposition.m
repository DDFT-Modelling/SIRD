function PG = Box_Decomposition(g, p,lb,ub, type)
% Decomposes vector g in three functions according to p:
% - Active:   p satisfies box contraints
%        L:   for the lower bound
%        U:   for the upper bound
% - Inactive: p does not satisfy box constraints

if nargin <= 4
  type = 0;
end

if ~isnan(type)
    AL = zeros(size(p)); %NaN if we want to display just
    AU = zeros(size(p));
    IN = zeros(size(p));
else
    AL = NaN(size(p)); %NaN if we want to display just
    AU = NaN(size(p));
    IN = NaN(size(p));
end
    
    % Active lower bound
    AL(p==lb) = g(p==lb);
    AU(p==ub) = g(p==ub);
    
    % Inactive components
    IN((p > lb) & (p< ub)) = g((p > lb) & (p< ub));
    
    % Return result
    %PG = [AL,AU,IN];
    PG(:,:,1) = AL;
    PG(:,:,2) = AU;
    PG(:,:,3) = IN;
end