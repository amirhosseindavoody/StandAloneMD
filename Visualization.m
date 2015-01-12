%%
% close all;
clear all;
clc;
fig = 30;
eV = 1.6e-19;
numAtom = 2;
directory = 'C:\Users\amirhossein\Documents\GitHub\StandAloneMD\StandAloneMD\bin\Release\';
%% Plot potential
% filename = 'potential.txt';
% myPotential = load ([directory filename]);
% [nR,m] = size(myPotential);
% 
% fig=fig+1; figure(fig);
% for i = 1:m
%     plot(myPotential(:,i),'-r','LineWidth',3); hold on;
% end;
% axis tight;
% return;

%% Plot force
% filename = 'force.txt';
% myForce = load ([directory filename]);
% [nR,m] = size(myForce);
% 
% fig=fig+1; figure(fig);
% for i = 1:m
%     plot(myForce(:,i),'-r','LineWidth',3); hold on;
% end;
% axis tight;
% return;

%% Plot pair distribution
filename = 'pairDistribution.txt';
pairDistribution = load ([directory filename]);
dR=0.1;
r = linspace(dR,length(pairDistribution)*dR,length(pairDistribution));

fig=fig+1; figure(fig);
plot(r,pairDistribution(:,1),'-b','LineWidth',3); hold on;
axis tight;
% return;

%% Plot energy
filename = 'energy.txt';
energy = load ([directory filename]);

fig=fig+1; figure(fig);
plot(energy(:,1)/eV,'-b','LineWidth',5); hold on;
plot(energy(:,2)/eV,'-r','LineWidth',5);
plot((energy(:,1)+energy(:,2))/eV,'-k','LineWidth',5);
axis tight;
return;

%% Plot temperature
% filename = 'temperature.txt';
% temperature = load ([directory filename]);
% 
% fig=fig+1; figure(fig);
% plot(temperature(:,1),'-b','LineWidth',3);
% axis tight;
% return;

%% visulize atom movements
filename = 'position.txt';
position = load ([directory filename]);
nTime = size(position,1)/numAtom;

Max = +10;
Min = -10;

fig = fig+1; figure(fig);
for iT = 1 :5: nTime
% for iT = nTime
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
%     pause (0.1);
    pause;
end;