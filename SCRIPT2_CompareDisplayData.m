% SCRIPT_DisplayData
%   Plot live optitrack data alongside rostopic data
%   Can be used for comparing 
%
%   TODO, add in mavros conversion
%   Some code adapted from https://github.com/kutzer/OptiTrackToolbox

%% Clear workspace, close all figures, clear command window
clear all
close all
clc
rosshutdown;

%% Create OptiTrack object and initialize
obj = OptiTrack;
Initialize(obj,'192.168.1.5','multicast');

%% Create Live Graph
figure;
%xlim([0 25]);
%ylim([-1.1 1.1]);
%x = 1:0.01:25;
%y = sin(x);

%% Subscribe to /uav0/mavros/setpoint_position/local to compare
rosinit('192.168.1.12');
% Ensure that subscribed topic is broadcasting else script will freeze
setpointsub=rossubscriber('/uav0/mavros/vision_pose/pose',"DataFormat","struct");
%setpointsub=rossubscriber('/uav0/mavros/setpoint_position/local',"DataFormat","struct");

%% Variables
bodyname='VisionVroom';

px=0; % previous location of body, used to plot lines
py=0; 
pz=0;

px2=0;
py2=0;
pz2=0;

init=0;

%% Display data
while true

    rb = obj.RigidBody; % Get current Optitrack rigid body information
    posedata = receive(setpointsub,10); % Get ROS pose information

    % Output frame information
    fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);
    
    % Plot origin
    plot(0,0,'+k');
    
    % Update each rigid body
    for i = 1:numel(rb)
        
        % Check for correct body
        % This crashes if there is a glitching body, remove any unused body from the asset tab
        if rb(i).Name == bodyname 
        
            % If body loses tracking, stop plotting
            if isempty(rb(i).Position)== 0
                
                if init ~= 0
                    fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
                    fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                    %plot(rb(i).Position(1)/1000,rb(i).Position(2)/1000,'-or');
                    plot3([rb(i).Position(1)/1000 px/1000],[rb(i).Position(2)/1000 py/1000],[rb(i).Position(3)/1000 pz/1000],'b');
                    plot3([posedata.pose.position.x px2], [posedata.pose.position.y py2], [posedata.pose.position.z pz2],'r');
                    
                    px = rb(i).Position(1);
                    py = rb(i).Position(2);
                    pz = rb(i).Position(3);
                    
                    px2 = posedata.pose.position.x;
                    py2 = posedata.pose.position.y;
                    pz2 = posedata.pose.position.z;
                    
                % On initalisation, there is no previous point to plot a line towards
                else
                   disp("initalising new plot")
                   
                    px = rb(i).Position(1);
                    py = rb(i).Position(2);
                    pz = rb(i).Position(3);
                    
                    px2 = posedata.pose.position.x;
                    py2 = posedata.pose.position.y;
                    pz2 = posedata.pose.position.z;
                   
                    plot3([rb(i).Position(1)/1000 px/1000],[rb(i).Position(2)/1000 py/1000],[rb(i).Position(3)/1000 pz/1000], 'b');
                    plot3([posedata.pose.position.x px2], [posedata.pose.position.y py2], [posedata.pose.position.z pz2],'r');
                   init=1;
                end
                %title('Title')
                xlabel('x') 
                ylabel('y')
                zlabel('z') 
                grid on;
                hold on;
                drawnow;
            end
        end
    end
end
        
        
        