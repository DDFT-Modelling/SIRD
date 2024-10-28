%%
% First contour plot in [0,1]^2
msize = 60;
F = zeros(msize,msize);
C = linspace(0,1,msize);
W = linspace(0,1,msize);

for i = 1:msize
    for j = 1:msize
        State_t = State(C(i),W(j), rho_ic, dims, box, aLine,tols);
        f = Objective(State_t, Target_t, N, Int_Spatial, Int_Time)/(N_T^2);
%         if W(j) == 0
%             display([W(j),f])
%         end
%         if W(j) == 1
%             display([W(j),f])
%         end
        F(i,j) = f;
    end
end

[X,Y] = meshgrid(C,W);

%[Cont,h] = contourf(X,Y,F',20);
%[Cont,h] = contourf(X,Y,F',[1, 0.7,0.4,0.3, 0.2, 5e-2, 5e-3, 1e-3, 0]);
[Cont,h] = contourf(X,Y,F',[0.25, 0.2, 0.14, 0.1, 0.06, 0.03, 0.015, 0.004 0.001, 0.0002, 0]);
h.EdgeColor = [142,140,169]/255;
xlabel('Contact rate $\beta$','Interpreter','latex','fontsize',11,'color','#555555')
ylabel('Recovery rate $\gamma$','Interpreter','latex','fontsize',11,'color','#555555')
set(gca,'XColor','#555555')
set(gca,'YColor','#555555')
colormap(viridis)

cb = colorbar();
cb.Ruler.Exponent = -2;
cb.Color = '#555555';
cb.Box = 'off';
cb.TickDirection = 'out';
cb.TickLabelInterpreter = 'latex';


set(gca,'TickDir','out','Box','off');
box off
hold on
plot(0.479,0.125,'w.')
hold off

original_fig = figure(3);
original_fig.Position = [1073         484         444         280];
set(gca, 'FontSize', 12);

exportgraphics(figure(3), 'Diffusive_Contours_Coarse.pdf', 'BackgroundColor','none','ContentType','vector')


%%
% Second contour plot in [0.4,0.6] x [0.05,0.3]

msize = 40;
F = zeros(msize,msize);
C = linspace(0.4,0.6,msize);  % 0.479,0.125
W = linspace(0.05,0.3,msize);  %

for i = 1:msize
    for j = 1:msize
        State_t = State(C(i),W(j), rho_ic, dims, box, aLine,tols);
        f = Objective(State_t, Target_t, N, Int_Spatial, Int_Time)/(N_T^2);
        F(i,j) = f;
    end
end

figure(4)
[X,Y] = meshgrid(C,W);

%[Cont,h] = contourf(X,Y,F',10);
%[Cont,h] = contourf(X,Y,F',[0.04, 0.03, 0.015, 0.004 0.002 0.001, 0.0002, 0]);
[Cont,h] = contourf(X,Y,F',[0.25, 0.2, 0.14, 0.1, 0.06, 0.03, 0.015, 0.004 0.001, 0.0002, 0]);
h.EdgeColor = [142,140,169]/255;
xlabel('Contact rate $\beta$','Interpreter','latex','fontsize',11,'color','#555555')
ylabel('Recovery rate $\gamma$','Interpreter','latex','fontsize',11,'color','#555555')
set(gca,'XColor','#555555')
set(gca,'YColor','#555555')
colormap(viridis)

cb = colorbar();
cb.Ruler.Exponent = -2;
cb.Color = '#555555';
cb.Box = 'off';
cb.TickDirection = 'out';
cb.TickLabelInterpreter = 'latex';


set(gca,'TickDir','out','Box','off');
box off

hold on
plot(0.479,0.125,'w.')
hold off

original_fig = figure(4);
original_fig.Position = [626         484         444         280];
set(gca, 'FontSize', 12);

exportgraphics(figure(4), 'Diffusive_Contours_Finer.pdf', 'BackgroundColor','none','ContentType','vector')
