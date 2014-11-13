%%
close all;
clear all;
clc;
fig = 0;
eV = 1.6e-19;

%%
% filename = 'C:\Amirhossein\StandAloneMD\StandAloneMD\bin\Debug\force.txt';
% rawForce = load (filename);
% 
% nR = size(rawForce,1);
% 
% fig = fig+1; figure(fig);
% plot(rawForce,'LineWidth',3);

%%
filename = 'C:\Amirhossein\StandAloneMD\StandAloneMD\bin\Debug\position.txt';
position = load (filename);

numAtom = 100;
nTime = size(position,1)/numAtom;

Max = +5;
Min = -5;

%plot kinetic, potential, and total energy of system
fig=fig+1; figure(fig);
plot(position(:,4)/eV,'-b','LineWidth',3); hold on;
plot(position(:,5)/eV,'-r','LineWidth',3);
plot((position(:,4)+position(:,5))/eV,'-k','LineWidth',3);
axis tight;
return;
%plot temperature of the gas
fig=fig+1; figure(fig);
plot(position(:,6),'-b','LineWidth',3);
axis tight;
return;
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
    pause (0.01);
%     pause;
end;