function g = Gradient(State_t, Target_t, Adjoint_t, N, Int_Spatial, Int_Time)
% Masks
maskS = 1:N;
maskI = N+1:2*N;
maskR = 2*N+1:3*N;
% Retrieve involved components from state
S    = State_t(:,maskS);       I = State_t(:,maskI);
% Retrieve each component of the target
Ta_S = Target_t(:,maskS);   Ta_I = Target_t(:,maskI);   Ta_R = Target_t(:,maskR);
% Retrieve each component of the adjoint
q_S  = Adjoint_t(:,maskS);   q_I = Adjoint_t(:,maskI);   q_R = Adjoint_t(:,maskR);

% Partial derivatives
dbeta  = Int_Time * (Int_Spatial * (S .* I .* ( q_I - q_S ))')';
dgamma = Int_Time * (Int_Spatial * (I .* ( q_R - q_I ))')';
% Gradient
g = [dbeta,dgamma];
end