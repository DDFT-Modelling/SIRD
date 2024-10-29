figure(500)

% First interpolation
IP = aLine.ComputeInterpolationMatrixPhys(0).InterPol;
rho_t = (IP * State_p)'/N_T;


Tile = tiledlayout(3,5,'TileSpacing','compact');

% Axes
for i = 1:5
h_S(i) = nexttile;
end
for i = 1:5
h_I(i) = nexttile;
end
for i = 1:5
h_R(i) = nexttile;
end

% Initial plot
axes(h_S(1))
box.plot(rho_t(maskS),{});
xlabel(''); ylabel('$x_2$'); caxis([0,2.5e-2])
title('$t=0$','Interpreter','latex');    shading interp

axes(h_I(1))
box.plot(rho_t(maskI),{});
xlabel(''); ylabel('$x_2$'); caxis([0,2.5e-2]);    shading interp

axes(h_R(1))
box.plot(rho_t(maskR),{});
xlabel('$x_1$'); ylabel('$x_2$'); caxis([0,2.5e-2]);    shading interp

set(h_S(1),'fontsize',6,'linewidth',0.1,'xticklabel',[]);
set(h_I(1),'fontsize',6,'linewidth',0.1,'xticklabel',[]);
set(h_R(1),'fontsize',6,'linewidth',0.1);

h_S(1).CameraPosition = [0 0 1];
h_I(1).CameraPosition = [0 0 1];
h_R(1).CameraPosition = [0 0 1];

% Second-Fourth interpolations
TimesEval = [0,3,5,15,30];

for i = 2:5
    t = TimesEval(i);
    IP = aLine.ComputeInterpolationMatrixPhys(t).InterPol;
    rho_t = (IP * State_p)'/N_T;
    
    % Add panels
    axes(h_S(i))
    box.plot(rho_t(maskS),{});
    xlabel(''); ylabel('');
    zlim([-1,1]); caxis([0,2.5e-2]);    shading interp
    title(['$t$ = ' num2str(TimesEval(i))],'Interpreter','latex')
    
    axes(h_I(i))
    box.plot(rho_t(maskI),{});
    xlabel(''); ylabel('');
    zlim([-1,1]); caxis([0,2.5e-2]);    shading interp

    axes(h_R(i))
    box.plot(rho_t(maskR),{});
    xlabel('$x_1$'); ylabel('');
    zlim([-1,1]); caxis([0,2.5e-2]);    shading interp;

    set(h_S(i),'fontsize',6,'linewidth',0.1,'xticklabel',[],'yticklabel',[]);
    set(h_I(i),'fontsize',6,'linewidth',0.1,'xticklabel',[],'yticklabel',[]);
    set(h_R(i),'fontsize',6,'linewidth',0.1,'yticklabel',[]);

    h_S(i).CameraPosition = [0 0 1];
    h_I(i).CameraPosition = [0 0 1];
    h_R(i).CameraPosition = [0 0 1];
end

% Add additional comments
set(h_S(5),'YAxisLocation', 'right')
set(h_I(5),'YAxisLocation', 'right')
set(h_R(5),'YAxisLocation', 'right')

axes(h_S(5));    ylabel('$\rho_S$');
axes(h_I(5));    ylabel('$\rho_I$');
axes(h_R(5));    ylabel('$\rho_R$');

colormap bone

exportgraphics(figure(500), strcat('Optimised_Densities.pdf'), 'BackgroundColor','none','ContentType','vector')

%colormap default

%cb = colorbar;
%cb.Layout.Tile = 'east';
%cb.Position = [0.92,0.12,0.02,0.8];


% Get measurements
original_fig = figure(500);
% original_position = original_fig.Position;        % 1056 420 560 318
original_fig.Position = [1056 420 560 320]