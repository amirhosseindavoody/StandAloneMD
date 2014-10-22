%%
close all;
clear all;
clc;
fig = 0;

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

numAtom = 10;
nTime = size(position,1)/numAtom;

Max = max(position,1);
Min = min(position,1);

fig = fig+1; figure(fig);
for iT = 1 : nTime
    
    iT
    for iAtom = 1:numAtom
        plot3(position((iT-1)*numAtom+iAtom,1),position((iT-1)*numAtom+iAtom,2),position((iT-1)*numAtom+iAtom,3),'b*'); hold on;
    end;
    xlim([0,10]);
    ylim([0,10]);
    zlim([0,10]);
    box on;
    hold off;
    pause;
end;