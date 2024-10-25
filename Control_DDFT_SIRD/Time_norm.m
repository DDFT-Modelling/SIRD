function a = Time_norm(Control, p, Int_Time)
    % ** Extract components **
    u = abs(Control(:,1));
    v = abs(Control(:,2));
    
    % ** Select p-norm **
    if length(p) == 1
        p = p * ones([2,1]);            % Scalar provided is converted into vector
    end
    p_1 = p(1);     % Inner vector
    p_2 = p(2);     % Time integral
    
    % ** Compute the norm **
    % Compute vector norm
    if p_1 < inf
        a_1 = ( u.^p_1 + v.^p_1 ).^(1.0/p_1);
    else
        a_1 = max(u,v);
    end
    % Compute integral norm in time
    if p_2 < inf
        a = (Int_Time * a_1.^p_2 ).^(1.0/p_2);
    else
        a = max(a_1,[],'all')';
    end
end