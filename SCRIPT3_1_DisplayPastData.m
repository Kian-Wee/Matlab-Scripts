% SCRIPT_DisplayPastData
% Plot historical data

% %% Clear workspace, close all figures, clear command window
% clear all
% close all
% clc

array_format = table2array(combinedposholdv2S7); % Change to Array Name of imported matrix

%% 1)Plot 2D/3D Data

% x deg, y deg, rps seperate
figure;
tiledlayout(3,1);

nexttile
plot(double(array_format(:,16)),double(array_format(:,5)),'linewidth',5);
% title('Body roll for position hold in clockwise direction with backmotor','fontweight','bold','fontsize',24);
title('Body roll for position hold in counter-clockwise direction with frontmotor','fontweight','bold','fontsize',24);
xlabel('Time','fontsize',24);
ylabel('Body roll in deg','fontsize',24);
grid on;

nexttile
plot(double(array_format(:,16)),double(array_format(:,6)),'linewidth',1,'color','r');
% title('Body pitch for position hold in clockwise direction with backmotor','fontweight','bold','fontsize',24);
title('Body pitch for position hold in counter-clockwise direction with frontmotor','fontweight','bold','fontsize',24);
xlabel('Time','fontsize',24);
ylabel('Body pitch in deg','fontsize',24);
grid on;

nexttile
plot(double(array_format(:,16)),double(array_format(:,19)),'linewidth',1,'color','k');
% title('Body RPS for position hold in clockwise direction with backmotor','fontweight','bold','fontsize',24);
title('Body RPS for position hold in counter-clockwise direction with frontmotor','fontweight','bold','fontsize',24);
xlabel('Time','fontsize',24);
ylabel('RPS','fontsize',24);
grid on;

hold off;