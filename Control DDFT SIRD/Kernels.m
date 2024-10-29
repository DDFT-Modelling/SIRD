function g = Kernels(y1,y2,data)
    % Compute the two kernels involved in the nonlocal interaction
    %
    % K_sd = exp( - \sigma_sd ( y1^2 + y2^2 ) )    % Social distancing
    % K_si = exp( - \sigma_si ( y1^2 + y2^2 ) )    % Self isolation
    %

    g1 = exp( - data.s_sd * (y1.^2 + y2.^2) );
    g2 = exp( - data.s_si * (y1.^2 + y2.^2) );
    g = [g1 g2];
end