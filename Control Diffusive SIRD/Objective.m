function f = Objective(State_t, Target_t, N, Int_Spatial, Int_Time)
% Difference of squares
Rho_t = (State_t - Target_t).^2;
% Integral
f = 0.5 * Int_Time * (Int_Spatial * ( Rho_t(:,1:N) + Rho_t(:,N+1:2*N) + Rho_t(:,2*N+1:end) )')';
end