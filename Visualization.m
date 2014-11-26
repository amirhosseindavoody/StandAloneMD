%%
% close all;
clear all;
clc;
fig = 20;
eV = 1.6e-19;
numAtom = 1000;
%%
% filename = 'C:\Amirhossein\StandAloneMD\StandAloneMD\bin\Debug\force.txt';
% rawForce = load (filename);
% 
% nR = size(rawForce,1);
% 
% fig = fig+1; figure(fig);
% plot(rawForce,'LineWidth',3);

%% Plot pair distribution
filename = 'C:\Amirhossein\StandAloneMD\StandAloneMD\bin\Release\pairDistribution.txt';
pairDistribution = load (filename);

%plot kinetic, potential, and total energy of system
fig=fig+1; figure(fig);
plot(pairDistribution(:,1),'-b','LineWidth',3); hold on;
axis tight;
return;

%% Plot energy
filename = 'C:\Amirhossein\StandAloneMD\StandAloneMD\bin\Release\energy.txt';
energy = load (filename);

%plot kinetic, potential, and total energy of system
fig=fig+1; figure(fig);
plot(energy(:,1)/eV,'-b','LineWidth',3); hold on;
plot(energy(:,2)/eV,'-r','LineWidth',3);
plot((energy(:,1)+energy(:,2))/eV,'-k','LineWidth',3);
axis tight;
% return;

%% Plot temperature
filename = 'C:\Amirhossein\StandAloneMD\StandAloneMD\bin\Release\temperature.txt';
temperature = load (filename);

%plot temperature of the gas
fig=fig+1; figure(fig);
plot(temperature(:,1),'-b','LineWidth',3);
axis tight;
return;

%% visulize atom movements
filename = 'C:\Amirhossein\StandAloneMD\StandAloneMD\bin\Release\position.txt';
position = load (filename);
nTime = size(position,1)/numAtom;

Max = +15;
Min = -15;

fig = fig+1; figure(fig);
for iT = 1 :20: nTime
    
    iT
    for iAtom = 1:numAtom
        plot3(position((iT-1)*numAtom+iAtom,1),position((iT-1)*numAtom+iAtom,2),position((iT-1)*numAtom+iAtom,3),'b*'); hold on;
    end;
    xlim([Min,Max]);
    ylim([Min,Max]);
    zlim([Min,Max]);
    box on;
    view([0 90]);
    hold off;
    pause (0.1);
%     pause;
end;