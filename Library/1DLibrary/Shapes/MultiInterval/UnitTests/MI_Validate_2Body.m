close all
clear all

yMin = -2;
yMax = 3;

opts1.endpoints = [yMin;yMax];
opts1.N = 50;
opts1.doPlots = false;

opts2.endpoints = [yMin;-1;1;yMax];
opts2.N = 50;
opts2.doPlots = false;

tic
[MI1,rho1] = MI_AdvectionDiffusion2Body(opts1);
time1 = toc;

tic
[MI2,rho2] = MI_AdvectionDiffusion2Body(opts2);
time2 = toc;

dy = 0.01;

interp1 = MI1.ComputeInterpolationMatricesUniform(dy);
interp2 = MI2.ComputeInterpolationMatricesUniform(dy);

nTimes = size(rho1,1);

err = 0;

pts1 = MI1.Pts.Uniform.y;
[pts2,keepMask,~] = unique(MI2.Pts.Uniform.y);

figure('Position',[100,100,1200,800])
for iTime = 1:nTimes
    rho1t = rho1(iTime,:)';
    rho2t = rho2(iTime,:)';
    rho1t = interp1*rho1t;
    rho2t = interp2*rho2t;
    rho2t = rho2t(keepMask);
    err = max(err,max(abs(rho1t-rho2t)));
    
    subplot(1,2,1)
    hold off
    plot(pts1,rho1t,'r');
    hold on
    plot(pts1,rho2t,'--b');
    subplot(1,2,2)
    plot(pts1,rho1t-rho2t,'r');
    pause(0.05)
end

disp(['Time 1 = ' num2str(time1)]);
disp(['Time 2 = ' num2str(time2)]);
disp(['Max error = ' num2str(err)]);