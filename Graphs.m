%% Clear workspace, close all figures, clear command window
%clear all
%close all
%clc
close all
% %% Plot
% 
% % plot3(t2.XPosition,t2.YPosition,t2.ZPosition,'b','LineWidth',3);
% plot3(t2.XSetpoint,t2.YSetpoint,t2.ZSetpoint,'b','LineWidth',3);
% % online GP + MPC
% % nominal MPC(first data)
% grid on
% title('Trajectory','fontname','times');
% xlabel("x(m)",'fontname','times');
% ylabel("y(m)",'fontname','times');
% zlabel("z(m)",'fontname','times');
% 
% figure;
% plot(t2.Timestamp,t2.EucaledianError,'LineWidth',3);
% grid on
% title('Position Error Magnitude against Time','fontname','times');
% ylabel("Position(m)",'fontname','times');
% xlabel("Time(s)",'fontname','times');
% 
% figure; hold on
% grid on
% a1 = plot(t2.Timestamp,t2.XVelocity,'LineWidth',3); M1 = "X-Velocity Error";
% a2 = plot(t2.Timestamp,t2.YVelocityFollowing,'LineWidth',3); M2 = "Y-Velocity Error";
% a3 = plot(t2.Timestamp,t2.ZVelocityFollowing,'LineWidth',3); M3 = "Z-Velocity Error";
% title('Velocity Error against Time','fontname','times');
% ylabel("Velocity(m/s)",'fontname','times');
% xlabel("Time(s)",'fontname','times');
% legend([a1, a2, a3], [M1, M2, M3]);
% hold off

%% Plot

% plot3(t2.XPosition,t2.YPosition,t2.ZPosition,'b','LineWidth',3);
plot3(m.XSetpoint,m.YSetpoint,m.ZSetpoint,'b','LineWidth',3);
% online GP + MPC
% nominal MPC(first data)
grid on
title('Trajectory','fontname','times');
xlabel("x(m)",'fontname','times');
ylabel("y(m)",'fontname','times');
zlabel("z(m)",'fontname','times');

figure('name','Moving Average'); hold on
grid on
a1 = plot(m.Timestamp,m.MovingAverage,':','LineWidth',4); M1 = "GPMPC";
a2 = plot(n.Timestamp,n.MovingAverage,'--','LineWidth',4); M2 = "Nominal MPC";
a3 = plot(pid.Timestamp,pid.MovingAverage,'-','LineWidth',4); M3 = "PID";
title('Euclidean Position Error against Time','fontname','times','fontsize',25);
ylabel("Error(m)",'fontname','times','fontsize',25, 'FontWeight', 'bold');
xlabel("Time(s)",'fontname','times','fontsize',25, 'FontWeight', 'bold');
legend([a1, a2, a3], [M1, M2, M3],'fontsize',25);
hold off

figure('name','Eucaledian Error'); hold on
grid on
b1 = plot(m.Timestamp,m.EucaledianError,'LineWidth',3); M11 = "GPMPC";
b2 = plot(n.Timestamp,n.EucaledianError,'LineWidth',3); M22 = "Nominal MPC";
b3 = plot(pid.Timestamp,pid.EucaledianError,'LineWidth',3); M33 = "PID";
title('Position Magnitude Error against Time','fontname','times');
ylabel("Position(m)",'fontname','times');
xlabel("Time(s)",'fontname','times');
legend([b1, b2, b3], [M11, M22, M33]);
hold off

figure('name','Velocity Error'); hold on
grid on
a1 = plot(m.Timestamp,m.XVelocityError,'-o','LineWidth',1); M1 = "X-Velocity Error";
a2 = plot(m.Timestamp,m.YVelocityError,'LineWidth',1); M2 = "Y-Velocity Error";
a3 = plot(m.Timestamp,m.ZVelocityError,'--','LineWidth',1); M3 = "Z-Velocity Error";
title('Velocity Error against Time','fontname','times');
ylabel("Velocity(m/s)",'fontname','times');
xlabel("Time(s)",'fontname','times');
legend([a1, a2, a3], [M1, M2, M3]);
hold off

%%
%figure;
%plot(t2.Timestamp,[t2.XSetpoint,t2.XPosition]);
%title('X Position with time');
%ylabel("Position(m)");
%xlabel("Time(s)");
%stackedplot(t2.XSetpoint,t2.XPosition);
%stackedplot(t2);

% figure;
% plot(t2.Timestamp,[t2.XVelocity,t2.XVelocityFollowing]);
% title('X-Velocity Error against Time','fontname','times');
% ylabel("Velocity(m/s)",'fontname','times');
% xlabel("Time(s)",'fontname','times');

% figure; hold on
% grid on
% a1 = plot(t2.Timestamp,t2.XVelocity); M1 = "X Chaser Velocity";
% a2 = plot(t2.Timestamp,t2.XVelocityFollowing); M2 = "X Target Velocity";
% title('X-Velocity Error against Time','fontname','times');
% ylabel("Velocity(m/s)",'fontname','times');
% xlabel("Time(s)",'fontname','times');
% legend([a1,a2], [M1, M2]);

% plot(t2.Timestamp,t2.XVelocityError);
% grid on
% title('X-Velocity Error against Time','fontname','times');
% ylabel("Velocity(m/s)",'fontname','times');
% xlabel("Time(s)",'fontname','times');
% 
% figure;
% plot(t2.Timestamp,t2.YVelocityError);
% grid on
% title('Y-Velocity Error against Time','fontname','times');
% ylabel("Velocity(m/s)",'fontname','times');
% xlabel("Time(s)",'fontname','times');
% 
% figure;
% plot(t2.Timestamp,t2.ZVelocityError);
% grid on
% title('Z-Velocity Error against Time','fontname','times');
% ylabel("Velocity(m/s)",'fontname','times');
% xlabel("Time(s)",'fontname','times');

% fig=figure  
% destination='C:\Users\MSI User\Desktop';
% print([destination,num2str(k),'.dpng']);
% close(fig)

% FolderName = 'C:\Users\MSI User\Desktop';   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   savefig(fullfile(FolderName, [FigName '.fig']));
% %   saveas(gcf,'Barchart.png');
%   export_fig(sprintf('figure%d',i),'-jpg');
% end

%%
% Nominal mpc
% Online GPMPC
% Plot all 3 trajectories