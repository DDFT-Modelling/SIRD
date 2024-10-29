function f = Objective(u, v, State_t, tik, N, Int_Spatial, Int_Time)
% Difference of squares
I = State_t(:,N+1:2*N);
% Integral
f = 0.5 * Int_Time * ((Int_Spatial * (I.^2)' )' + tik * (u.^2 + v.^2));
end

