% 
% 
% 
% 
% 
% 
% —————————————————————————————————————————————————————————————————————————
% 
% 
%  ¡¡¡¡The code in the inner loop should be updated from the notebook!!!!
%      Run up to line 104.
% 
% 
% —————————————————————————————————————————————————————————————————————————
% 
% 
% 
% 
% 
% 
% 


% 10 experiments, 2 curves
P = zeros(aLine.N,10,2);
pd_I = zeros(aLine.N,10);






Facts = [[10], 9:-1:1];
%Facts = [10,8,6,4,2,1];

%% Optimisation process
%------------------------------------------
%------------------------------------------
for index = 1:10


fact    = Facts(index);
Box_Lim = fact * D;

% Determine bounds and define projection
lb   = -Box_Lim;
ub   = min(Box_Lim, 10*D);                          %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Proj = @(x) min( max(x,lb), ub);


%Proj = @(x) min( max(x,-Box_Lim), Box_Lim);
%Now we can write a code for the algorithm:
Nb_Its = 100;       % Number of iterations
critical = 1;      % Report status after this number of iterations

% Define majorant function
quad  =@(x,y, h_y, dh_y, L) h_y + sum((x-y) .* dh_y,'all') + 0.5 * L * sumsqr(x-y);
quadS =@(x,y, h_y, dh_y, L) h_y + sum(Int_Time * ((x-y) .* dh_y)) + ...
                            0.5 * L * Time_norm(x-y, 2, Int_Time)^2;

% Define objects to store information across iterates
crit    = zeros(Nb_Its + 1, 1);
State_t = State_Init;
alpha   = [u_Init, v_Init];
%State_t = State_p;
%alpha   = p;

f       = Objective(u_Init, v_Init, State_t, tik, N, Int_Spatial, Int_Time);
f       = f/(N_T^2); % We scale the system
crit(1) = f;
time    = zeros(Nb_Its,1);      % Time per iteration
dnIn    = zeros(Nb_Its, 1);     % Relative norm of inactive gradient
rnCo    = zeros(Nb_Its, 1);     % Relative sq-norm of control
rdCo    = zeros(Nb_Its, 1);     % Relative norm of difference between controls
cdiff   = zeros(Nb_Its, 1);     % Absolute value two consecutive objective iterations
U_Store = zeros(n_t,Nb_Its+1);  % Iterations in u
V_Store = zeros(n_t,Nb_Its+1);  % Iterations in v
U_Store(:,1) = alpha(:,1);
V_Store(:,1) = alpha(:,2);

% Stopping criteria
Stop_norm = 1e-10;
Stop_crit = 1e-10;

fprintf('   k |   f(αₖ)  | ᶦ‖df(ωₖ)‖ᵣ|   Time   |   ‖αₖ‖ᵣ² | ‖αₖ-αₖ₋₁‖ᵣ| |fₖ-fₖ₋₁|  \n')
disp(strcat( repmat('–',1,75) ))

% Initialisation
omega = alpha;
theta = 1;      % Inertia parameter
eta = 1.1;      % Scaling for backtracking                    % was 10
L = 1.1^6;      % Approximation of the Lipschitz constant     % was 10
I = 50;         % Number of iterations for backtracking

tols = 1e-7;
timing = 0.0;
L_tol = 2e-3;
reduce_L = true;    % Activate reduction of estimate for Lipschitz constant
reducedL = false;   % Check if constant has been reduced already
% Iterations
for it = 1:Nb_Its
    % Start measuring time
    start = tic;
    
    alpha_old = alpha;
    t_old = theta;
    f_old = f;
    
    %---------------------------------------------------------------
    % Function evaluation and adjoint for omega
    State_o   = State(omega(:,1),omega(:,2), rho_ic, dims, aLine, Conv, Diff, tols);
    State_o   = max(State_o, 0.0); % non-negativity
    Adjoint_o = Adjoint(omega(:,1),omega(:,2), State_o, dims, aLine, Conv, Diff, tols);
        
    h_y  = Objective(omega(:,1),omega(:,2), State_o, tik, N, Int_Spatial, Int_Time)/(N_T^2);
    dh_y = Gradient(omega(:,1),omega(:,2), State_o, Adjoint_o, N, Int_Spatial, ...
        gConv_u, gConv_v, grad, G, tik)/(N_T^2);
    
    if norm(dh_y,'inf') > (Box_Lim * 1.1)
        scale = (Box_Lim * 1.1)/norm(dh_y);
        %scale = 1/norm(dh_y);
    else
        scale = 1.0;
    end
    
    if (it > 2) && reduce_L
        if cdiff(it) < L_tol    %abs(f_old - crit(it-1)) < L_tol
            L = max(min(0.1, L/10), 1e-3);
            L_tol = L_tol * 1e-1
        end
        reducedL = true;

%         if abs(f_old - crit(it-1)) < 1e-5
%             L = 0.01;
%         end
%         if abs(f_old - crit(it-1)) < 1e-3 * crit(it-1)
%             L = L;
%         end
    end

    %---------------------------------------------------------------
    % backtracking
    for i = 0:I
        % temporary approximation of Lipschitz constant
        L_tmp = eta^i * L;
        b = 1/L_tmp;
        % Gradient and proximal step
        q = omega - b * scale * dh_y;       % did i made a mess somewhere?
        %fprintf('u = (%.6f,%.6f) \n',u)
        p = Proj(q);
        % Compute state to test validation of L_tmp
        State_p = State(p(:,1),p(:,2), rho_ic, dims, aLine, Conv, Diff, tols);
        State_p = max(State_p, 0.0); % non-negativity
        h_p     = Objective(p(:,1),p(:,2), State_p, tik, N, ...
                                Int_Spatial, Int_Time)/(N_T^2);
        %fprintf('(h_y, h_p) = (%.6f,%.6f) \n',[h_y, h_p])
        %fprintf('(h, q) = (%.6f,%.6f) \n',[h_p,quadS(p,omega, h_y,dh_y, L_tmp/scale)])
        %fprintf('%.0f \n',[h_p<quadS(p,omega, h_y,dh_y, L_tmp/scale)])
        
        %if h_p <= quad(p,omega, h_y,dh_y, L_tmp) + (1e-2)
        if h_p <= quadS(p,omega, h_y,dh_y, L_tmp/scale) %+ (1e-2)
            break
        elseif reducedL
            reduce_L = false; % Prevent reduction to take place
        end
        fprintf('%d ',i)
        %fprintf('%d %f %f\n',i,h_p, quad(p,omega, h_y,dh_y, L_tmp))
    end
    %---------------------------------------------------------------
    % Update control
    L     = L_tmp;
    alpha = p;
    % inertia update
    theta = (it + 2.2)/(2.2);
    omega = alpha + ((t_old - 1)/theta) * (alpha - alpha_old);
    
    % Measure time up to this point!
    timed = toc(start);
    timing = timing + timed;
    % Compute derivative information
    PG = Box_Decomposition(dh_y, p, lb,ub); % Decompose gradient
    
    % Store information
    % --------------------------------------------------------------
    time(it+1)  = timing;
    crit(it+1)  = h_p;
    dnIn(it+1)  = Time_norm(PG(:,:,3), inf, Int_Time)/sqrt(nnz(PG(:,:,3)) );
    rnCo(it+1)  = Time_norm(alpha, 2, Int_Time)^2 / n_t;
    rdCo(it+1)  = Time_norm(alpha-alpha_old, 2, Int_Time) / sqrt(n_t);
    cdiff(it+1) = abs(h_p-f_old);
    U_Store(:,it+1) = alpha(:,1);
    V_Store(:,it+1) = alpha(:,2);
    
    % Display if needed
    if(mod(it,critical)==0)
        Progress = sprintf('%4.0f | %.3e | %.3e | %8.2f | %.3e | %.3e | %.3e', ...
            [it, h_p, dnIn(it+1), timing, rnCo(it+1), rdCo(it+1), cdiff(it+1)]);
        display(regexprep( Progress, '(?<=e[-+])0', '' ))
    end
    % --------------------------------------------------------------
    
    if min(State_p,[],'all') < 0.0
        fprintf('Negative! %3.3f %', 100*numel(State_p(State_p < 0.0))/numel(State_p))
    end
    
    if rdCo(it+1) < Stop_norm
        fprintf('‖αₖ - αₖ₋₁‖ᵣ < %.e\n', Stop_norm)
        
        time  = time(1:it+1);   crit  = crit(1:it+1);   dnIn  = dnIn(2:it+1);
        rnCo  = rnCo(2:it+1);   rdCo  = rdCo(2:it+1);   cdiff = cdiff(2:it+1);
        U_Store = U_Store(:,1:it+1);
        V_Store = V_Store(:,1:it+1);
        break
    end
    
    if cdiff(it+1) < Stop_crit
        fprintf('|f(αₖ) - f(αₖ₋₁)| < %.e\n', Stop_crit)
        
        time  = time(1:it+1);   crit  = crit(1:it+1);   dnIn  = dnIn(2:it+1);
        rnCo  = rnCo(2:it+1);   rdCo  = rdCo(2:it+1);   cdiff = cdiff(2:it+1);
        U_Store = U_Store(:,1:it+1);
        V_Store = V_Store(:,1:it+1);
        break
    end
    
    f = h_p;
    
end
disp(strcat( repmat('–',1,75) ))
fprintf('\nNorm of inactive gradient ᶦ‖df(αₖ)‖ ≈ ᶦ‖df(ωₖ)‖ = %.3e.\n\n', dnIn(it) )
fprintf('\nAlgorithm stopped after %.4f seconds and %.0f iterations.\n\n', ...
    [time(it+1),it])

    



P(:,index,1) = p(:,1);
P(:,index,2) = p(:,2);
pd_I(:,index) = Int_Spatial * State_p(:,maskI)'/N_T;

end
%------------------------------------------
%------------------------------------------





%% Plot controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal controls

figure(200)
% Axes
h_u = subplot(1,2,1);    h_v = subplot(1,2,2);

indices = [1,3,5,7,9,10];
% Plot u
axes(h_u)
hCurves_u = plot(outTimes, P(:,indices,1),'linewidth', 1.5);
ylim([-10*D,10*D]);    set(h_u,'fontsize',10,'linewidth',0.1);
title('$u(t)$ - Social distancing','Interpreter','latex','fontsize',15)
xlabel('Time $t$','Interpreter','latex','fontsize',12);
ylabel('Strength','Interpreter','latex','fontsize',12);

% Plot v
hold on
axes(h_v)
hCurves_v = plot(outTimes, P(:,indices,2),'linewidth', 1.5);
ylim([-10*D,10*D]);    set(h_v,'fontsize',10,'linewidth',0.1,'yticklabel', []);
title('$v(t)$ - Self-isolation','Interpreter','latex','fontsize',15)
xlabel('Time $t$','Interpreter','latex','fontsize',12);

% Set transparency
for index = 1:6
    if index ~= 1
        hCurves_u(index).Color = [hCurves_u(index).Color, 0.3];
        hCurves_v(index).Color = [hCurves_v(index).Color, 0.3];
    end
end

% Set grids
grid(h_u,'on');    grid(h_v,'on');

% Legend
axes(h_u)
leg = legend(string( Facts(indices)*(-D) ),'Interpreter','latex','Location','NorthEast','fontsize',10);
title(leg,'$\ell_u$','Interpreter','latex','fontsize',12);
axes(h_v)
leg = legend(string( Facts(indices)*(-D) ),'Interpreter','latex','Location','NorthEast','fontsize',10);
title(leg,'$\ell_v$','Interpreter','latex','fontsize',12);

% position
h_u(1).Position = [ 0.13    0.11    0.36    0.75];
h_v(1).Position = [ 0.53    0.11    0.36    0.75];

sgtitle('Optimal interaction strengths','Interpreter','latex','FontSize',17)


exportgraphics(figure(200), strcat('Optimised_Controls_',num2str(tik),'.pdf'), 'BackgroundColor','none','ContentType','vector')
hold off




%% Plot aggregated curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effect in infected

% Select subvalues for zoom
sub_outTimes_Bool = (outTimes <= 25) & (outTimes >= 12);
subTimes = outTimes( sub_outTimes_Bool );
sub_pd_I = pd_I(sub_outTimes_Bool,: );
% Integrate infected without control
I_No_uv = Int_Spatial * State_Init(:,maskI)'/N_T;
sub_I_No_uv = I_No_uv(sub_outTimes_Bool);

% Start figure
figure(300)
A1 = axes();
set(A1,'fontsize',10,'linewidth',0.1);

hCurves_pd_I = plot(outTimes, pd_I(:,indices),'linewidth', 1.5);
hold on
plot(outTimes, I_No_uv, 'Color','black', 'LineStyle', '--')

xlabel('Time $t$','Interpreter','latex','fontsize',12);
ylabel('Infected population','Interpreter','latex','fontsize',12);
title('Evolution of infected by strategy','Interpreter','latex','fontsize',17)
hold off

for index = 1:6
    if index ~= 1
        hCurves_pd_I(index).Color = [hCurves_pd_I(index).Color, 0.3];
    end
end
% Legend
leg = legend(string( Facts(indices)*(-D) ),'Interpreter','latex','Location','NorthEast','fontsize',10);
title(leg,'$\ell$','Interpreter','latex','fontsize',12);

% Zoomed box
A2 = axes('position',[.3 .5 .4 .4 ]);
box on
sub_hCurves_pd_I = plot(subTimes, sub_pd_I(:,indices),'linewidth', 1.5);
hold on
plot(subTimes, sub_I_No_uv, 'Color','black', 'LineStyle', '--')

for index = 1:6
    if index ~= 1
        sub_hCurves_pd_I(index).Color = [sub_hCurves_pd_I(index).Color, 0.3];
    end
end

ylim([min(sub_pd_I,[],'all'), max(sub_pd_I,[],'all')])
xlim([subTimes(1), subTimes(end)])
set(A2,'fontsize',10,'linewidth',0.1);

% Activate grid lines
grid(A1,'on');    grid(A2,'on');

% Annotate zoom
annotation('rectangle',[0.324,0.169,0.19,0.13], 'Color', [0.5,0.5,0.5])
annotation('arrow',[0.4179 0.4179+0.042], [0.321 0.321+0.0986], 'Color', [0.5,0.5,0.5])

% Store plot
exportgraphics(figure(300), strcat('Optimised_Infected_',num2str(tik),'.pdf'), 'BackgroundColor','none','ContentType','vector')



%% Store data

save('store_Box_uv.mat','P');
save('store_Box_I.mat','pd_I');

%
P = matfile('store_Box_uv.mat').P;
pd_I = matfile('store_Box_I.mat').pd_I;





%% Plot densities
State_p   = State(P(:,1,1),P(:,1,2), rho_ic, dims, aLine, Conv, Diff, tols);
% Plot gif
plots_in_box(State_p,  N_T, maskS,maskI,maskR, box,outTimes,Int_Spatial)
% Plot aggregated compartments
Plot_SIR_Mean_Curves


% Set colours
% cols = 8;
% left_color = [102, 87, 139]/255;       % white
% right_color = [98, 215, 217]/255; % Soft blue
% cmap = interp1([0, 1], [left_color; right_color], linspace(0, 1, cols));
% 
% for index = 1:6
%     if index ~= 1
%         hCurves_u(index).Color = cmap(index,:);
%         hCurves_v(index).Color = cmap(index,:);
%     else
%         hCurves_u(index).Color = [0,0,0];
%         hCurves_v(index).Color = [0,0,0];
%     end
% end



















