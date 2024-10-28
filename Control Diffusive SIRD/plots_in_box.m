function plots_in_box(Rho_t, N_T, maskS,maskI,maskR, box,outTimes,Int_Spatial, filename)

if nargin <= 8
  filename = strcat( 'SIR-Forward.gif' );
end
    % Get masked result for each component
    S_t = Rho_t(:, maskS)/N_T;
    I_t = Rho_t(:, maskI)/N_T;
    R_t = Rho_t(:, maskR)/N_T;

    % constants for plotting
    %maxRho = max(max(Rho_t));
    %minRho = min(min(Rho_t));
    %m_ic = Int*rho_ic(N+1:2*N);

    % Plot solution
    h  = figure('Position',[100,100,1200,400]);
    axis tight manual % this ensures that getframe() returns a consistent size

    MaxS = max(S_t,[],'all');   MinS = min(S_t,[],'all');
    MaxI = max(I_t,[],'all');   MinI = min(I_t,[],'all');
    MaxR = max(R_t,[],'all');   MinR = min(R_t,[],'all');
    % Add a little extra of height with the same magnitude
    MaxS = MaxS + 10^floor( log10(abs(MaxS)) );
    MaxI = MaxI + 10^floor( log10(abs(MaxI)) );
    MaxR = MaxR + 10^floor( log10(abs(MaxR)) );

    hS = subplot(1,3,1);  zlim(hS,[MinS MaxS]);  caxis(hS,[0 MaxS]);
    hI = subplot(1,3,2);  zlim(hI,[MinI MaxI]);  caxis(hI,[0 MaxI]);
    hR = subplot(1,3,3);  zlim(hR,[MinR MaxR]);  caxis(hR,[0 MaxR]);

    set(gca,'Color','none');   set(gcf, 'color', 'white');

    opts = {};  % plotting options - the default are ok for now
    

    for iTime = 1:length(outTimes)

        % extract the solution at this time
        St = S_t(iTime,:)';
        It = I_t(iTime,:)';
        Rt = R_t(iTime,:)';
        I_int(iTime) = Int_Spatial*It;%sum(It*Int,'all');
        S_int(iTime) = Int_Spatial*St;
        R_int(iTime) = Int_Spatial*Rt;
        %rhot = rhot(n+1:2*n);  % If line 103 is commented

        % plot with built-in class function
        hold off
        % Plot S
        axes(hS)
        box.plot(St,opts);
        zlim(hS,[MinS MaxS]);   caxis(hS,[0 MaxS]);
        xlabel(''); ylabel(''); zlabel('')
        title('$S$','Interpreter','latex')
        hold off
        % Plot I
        axes(hI)
        box.plot(It,opts);
        zlim(hI,[MinI MaxI]);   caxis(hI,[0 MaxI]);
        xlabel(''); ylabel(''); zlabel('')
        title('$I$','Interpreter','latex')
        hold off
        % Plot R
        axes(hR)
        box.plot(Rt,opts);
        zlim(hR,[MinR MaxR]);   caxis(hR,[0 MaxR]);
        xlabel(''); ylabel(''); zlabel('')
        title('$R$','Interpreter','latex')


        % set title with some info
        sgtitle(['$t$ = ' num2str(outTimes(iTime))],'Interpreter','latex','FontSize',20)

        % fix axes
    %         colorbar

        %pause(0.05)

        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 

        % Write to the GIF File 
        if iTime == 1 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
        end 

    end
end