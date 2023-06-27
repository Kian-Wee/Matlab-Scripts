% SCRIPT1_2_StreamDataROS2
%   Stream live data received by from OptiTrack software(Motive) to ROS2

%% Clear workspace, close all figures, clear command window
clear all
close all
clc

import optitrack_class.*

%% Create OptiTrack object, ROS node and publisher
optitrack_node = ros2node("/matlab_optitrack_publisher");
obj = OptiTrack;
Initialize(obj,'192.168.1.5','multicast'); %Ensure broadcast frame id is on, loop interface is set to this ip and transmission type is set to multicast

%% Display data
prb = obj.RigidBody; % previous location of body, used to plot lines

desiredRate = 60; % approx 60 hz
rate = robotics.Rate(desiredRate);
rate.OverrunAction = 'drop';

% Initalise all bodies at first
rb = obj.RigidBody;
% Create each motive defined rigid body
arr = blanks(numel(rb));

for i = 1:numel(rb)
    my_field = strcat(convertCharsToStrings(rb(i).Name));
    variable.(my_field) = optitrack_class;
    variable.(my_field).init(convertCharsToStrings(rb(i).Name),optitrack_node,"px4");
end

disp('Outputing poses for: ');
disp(variable)

%% Loop

while true    
    % Get current rigid body information, this has to be recalled every time for new frame index
    rb = obj.RigidBody;
   
    % Output frame information
    % fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);
    % Update each rigid body
    for i = 1:numel(rb)
        % This crashes if there is a glitching body, remove any unused body
        % from the asset tab, if not uncomment this line
            if isempty(rb(i).Position)== 0 % Check that position array is properly formed(not missing)
                rb_position = rb(i).Position/1000;
                rb_quaternion = rb(i).Quaternion; % w, x, y, z
                   
%                 fprintf('%s\t   Position [%f,%f,%f]\n', rb(i).Name, rb_position);
%                 fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb_quaternion );
%                 hydra1_pose_msg.pose.position.x = rb_position(1);
%                 hydra1_pose_msg.pose.position.y = rb_position(2);            
%                 hydra1_pose_msg.pose.position.z = rb_position(3);
%                 hydra1_pose_msg.pose.orientation.w = rb_quaternion(1,1);
%                 hydra1_pose_msg.pose.orientation.x = rb_quaternion(1,2);
%                 hydra1_pose_msg.pose.orientation.y = rb_quaternion(1,3);
%                 hydra1_pose_msg.pose.orientation.z = rb_quaternion(1,4);

                my_field = strcat(convertCharsToStrings(rb(i).Name));
                variable.(my_field).sendpose(rb_position(1),rb_position(2),rb_position(3),rb_quaternion(1,1),rb_quaternion(1,2),rb_quaternion(1,3),rb_quaternion(1,4));
            end
%         else
%             arr(end+1) = rb(i).Name;
%             %create object   
    end
    waitfor(rate);
end