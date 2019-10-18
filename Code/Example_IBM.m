%% EXAMPLE IBM REALISATION

% Experiment indices to match raw data
repnames = {'8000_E2', '8000_G2', '10000_H1',  ...
            '10000_H2','8000_F2', '10000_A2', ...
            '12000_H3','12000_F3','12000_D2'};
        
% Experiment to load initial condition
ExpIndex = 5;

% Load Initial Condition
IC = csvread(['../Data/',repnames{ExpIndex},'/PC3_',repnames{ExpIndex},'_12h.csv']);

% Parameters
m  =  0.30;     % motility rate
p  =  0.03;     % proliferation rate
gm = -1.00;     % motility interaction strength
gp =  0.01;     % proliferation interaction strength
gb = 50.00;     % directional bias strength
s  = 12.00;     % sigma (interaction distance)
md = 24.00;     % movement step distance

params = [m,p,gm,gp,gb,s,md];

% Domain
L  = 1440;
H  = 1900;

% End time
T  = 36;

% Max number of agents in simulation
Nmax = 5000;

% Seed
seed = 1;

% Run model
[XT,X18] = IBM(params,[L,H],IC,T,Nmax,seed);

% Plot
clf;
Locations = {IC,X18,XT};
for i = 1:3
    
    % Plot Experimental Data
    subplot(3,2,i*2-1);
    Iexp = Plot_ScatterImage(repnames{ExpIndex},(i-1)*18+12);
    imshow(Iexp);
    title(['t = ',num2str((i-1)*18),' h']);
    
    % Plot Simulation Data
    subplot(3,2,i*2);
    Isim = Plot_ScatterImage(Locations{i});
    imshow(Isim);
    title(['t = ',num2str((i-1)*18),' h']);
    
end
