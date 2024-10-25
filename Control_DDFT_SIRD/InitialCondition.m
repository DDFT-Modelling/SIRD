function [y,N_T] = InitialCondition(y1,y2,Int,N)
choice = 'Gauss&Donut';     % Name of initial condition
f = 0.5;                    % Scale factor for infected

if choice == 'Gauss&Donut'
    % S initial condition is a gaussian with amplitude 7.964, variance 100/50 = 2, centred at 0
    % I initial is a mollified donut with volume (f * Int * S)
    % R initial is 0

    S = 7.964 * exp( -0.75 * ((y1-2).^2 + (y2+2).^2) );
    S = 0.7 * S + 7.964 * exp( -0.75 * ((y1+0.75).^2 + (y2-1.5).^2) ) * 0.3;
    S_T = Int * S; % Total number of people from susceptible

    I = exp( -1 *( (y1+0.75).^2 + (y2-1.5).^2 ));
    II = exp( -3 *( sqrt( (y1+2.5).^2 + (y2+3).^2) - 1).^2);
    I = 0.75*I + 0.25* II;
    I = S_T * I/(Int*I);

    % Scale
    I = f * I;
    S = (1-f) * S;

    R = zeros([N,1]);

elseif choice == 'TwoGauss'
    S = 7.964 * exp( -0.75 * (y1.^2 + y2.^2) );
    S_T = Int * S; % Total number of people from susceptible

    I = exp( -1 *( (y2-1.5).^2 + (y1+0.75).^2));
    I = S_T * I/(Int*I);
    
    I = f * I;
    S = (1-f) * S;

    R = zeros([N,1]);
    
elseif choice == 'Donut'
    S = 7.964 * exp( -0.75 * (y1.^2 + y2.^2) );
    S_T = Int * S; % Total number of people from susceptible

    I = exp( -1 *( sqrt(y2.^2 + y1.^2) - 5).^2); % 5
    I = S_T * I/(Int*I);

    I = f * I;          % for the other I: exp( -1 *( (y2-1.5).^2 + (y1+0.75).^2));
    S = (1-f) * S;

    R = zeros([N,1]);
    
elseif choice == 'Scaled'
	S = 7.964 * exp( -0.25 * (y1.^2 + y2.^2) );
    I = f * S;
    S = (1-f) * S;
    R = zeros([N,1]);  
end

N_T = sum([Int*S,Int*I,Int*R]);    % Compute total number of people in the population
y   = [S;I;R];
end