% Plot historical rosbag data

%% Clear workspace, close all figures, clear command window
clear all
close all
clc

%% Initalize variables
%filename='C:\%userprofile%\Desktop\1.bag'; %doesnt work
filename='C:\Users\MSI User\Desktop\1.bag';
%topic='/uav0/mavros/setpoint_position/local'; %Position of drone from mavros
topic='/uav0/mavros/vision_pose/pose'; %Position of drone from external vision system(eg optitrack)

%% Open bag file 
% Topic, DataFormat and struct are generic, change filename and topic
bag = rosbag(filename);
pos = select(bag, 'Topic', topic);
msgStructs = readMessages(pos,'DataFormat','struct');

%start = bag.StartTime
%bagselect2 = select(bag,"Time",[start start + 30],"Topic","/odom") %To read first 30 seconds
%[start, start+1; start+5, start+6]% to stich 2 timings
%disp(msgStructs{1}.Pose.Position)

%% 1)Plot 2D/3D Data
for i= 2:length(msgStructs)
    % 2D
    %plot(rb(i).Position(1)/1000,rb(i).Position(2)/1000,'-or');
    % 3D
    plot3([msgStructs{i-1}.Pose.Position.X msgStructs{i}.Pose.Position.X],[msgStructs{i-1}.Pose.Position.Y msgStructs{i}.Pose.Position.Y],[msgStructs{i-1}.Pose.Position.Z msgStructs{i}.Pose.Position.Z],'b');
end

xlabel('x');
ylabel('y');
zlabel('z');
title('3D Plot of object trajectory');
grid on;

%% 2)Plot 1D Data Against Time
% Where the first value is the array for all the values and the subsequent values are the rows(topics) to be displayed
ts=timeseries(pos,"Pose.Position.X","Pose.Position.Y","Pose.Position.Z");
figure
plot(ts,"LineWidth",3)
colorbar