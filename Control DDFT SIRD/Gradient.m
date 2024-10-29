function g = Gradient(u,v, State_t, Adjoint_t, N, Int_Spatial, gConv_u, gConv_v, grad, G, tik)
% Masks
maskS = 1:N;
maskI = N+1:2*N;
maskR = 2*N+1:3*N;
% Paramters
G_S = G{1};
G_I = G{2};
G_R = G{3};
% Retrieve involved components from state
S    = State_t(:,maskS);       I = State_t(:,maskI);       R = State_t(:,maskR);
% Retrieve each component of the adjoint
q_S  = Adjoint_t(:,maskS);   q_I = Adjoint_t(:,maskI);   q_R = Adjoint_t(:,maskR);


% Precompute vectorial parts
RHS_aux = [S,S]' .* (grad * (G_S * q_S)') + [R,R]' .* (grad * (G_R * q_R)');
RHSdu = RHS_aux .*  (gConv_u * (S + R)');
RHSdv = RHS_aux .*  (gConv_v * I')  +  ...
            ( [I,I]' .* (grad * (G_I * q_I)') )  .*  (gConv_v * (S + I + R)');

% Partial derivatives
du = tik * u + (Int_Spatial * (RHSdu(1:N,:) + RHSdu(N+1:2*N,:)))';
dv = tik * v + (Int_Spatial * (RHSdv(1:N,:) + RHSdv(N+1:2*N,:)))';

g = [du,dv];
end