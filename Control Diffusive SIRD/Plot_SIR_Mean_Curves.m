
S_int = Int_Spatial * State_p(:,maskS)';
I_int = Int_Spatial * State_p(:,maskI)';
R_int = Int_Spatial * State_p(:,maskR)';

[argvalue, argmax] = max(I_int/(max(S_int) + I_int(1) ));

hold on
plot(outTimes,S_int/(max(S_int) + I_int(1) ), 'Color',[0, 0.4470, 0.7410])
plot(outTimes,I_int/(max(S_int) + I_int(1) ), 'Color',[0.9290, 0.6940, 0.1250])
plot(outTimes,R_int/(max(S_int) + I_int(1) ), 'Color',[0.4940, 0.1840, 0.5560])

% First maximum
line([outTimes(argmax) outTimes(argmax)], [0 argvalue], 'Color','red','LineStyle','--');
% Second maximum
[argvalue, argmax] = max(I_int(30:end)/(max(S_int) + I_int(1) ));
line([outTimes(argmax+30) outTimes(argmax+30)], [0 argvalue], 'Color','red','LineStyle','--');

hl = legend('$\overline{S}(t)$','$\overline{I}(t)$','$\overline{R}(t)$');
set(hl, 'Interpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
xlabel('$t$','interpreter','latex')
ylabel('Total population fraction','interpreter','latex')
%title('Population fractions for each compartment','Interpreter','latex')
grid(gca,'on');

exportgraphics(figure(3), strcat('Integrated_Curves.pdf'), 'BackgroundColor','none','ContentType','vector')






