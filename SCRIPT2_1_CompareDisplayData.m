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
Initialize(obj,'192.168.1.25','multicast');

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
vx=0;
vy=0;
vz=0;

px2=0;
py2=0;
pz2=0;
vx2=0;
vy2=0;
vz2=0;

init=0;
ptime=0; % For velocity caculations
x=[]; % Array to store to excel

%% Break loop if keypress to save to excel
DlgH = figure;
H = uicontrol('Style', 'PushButton', ...
                    'String', 'Break', ...
                    'Callback', 'delete(gcbf)');

%% Display data
while (ishandle(H))

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
                    
                    vx = (rb(i).Position(1)-px)/(cputime-ptime)/1000;
                    vy = (rb(i).Position(2)-py)/(cputime-ptime)/1000;
                    vz = (rb(i).Position(3)-pz)/(cputime-ptime)/1000;
                    euv = sqrt(vx^2+vy^2+vz^2);
                    
                    vx2 = (posedata.pose.position.x-px)/(cputime-ptime)/1000;
                    vy2 = (posedata.pose.position.y-py)/(cputime-ptime)/1000;
                    vz2 = (posedata.pose.position.z-pz)/(cputime-ptime)/1000;
                    euv2 = sqrt(vx2^2+vy2^2+vz2^2);
                    
                    px = rb(i).Position(1);
                    py = rb(i).Position(2);
                    pz = rb(i).Position(3);
                    
                    px2 = posedata.pose.position.x;
                    py2 = posedata.pose.position.y;
                    pz2 = posedata.pose.position.z;
                    
                    ptime=cputime;
                    
                    x=[x; [rb(i).Position(1) rb(i).Position(2) rb(i).Position(3) vx vy vz euv
                        posedata.pose.position.x posedata.pose.position.y posedata.pose.position.z vx2 vy2 vz2 euv2]];
                    
                % On initalisation, there is no previous point to plot a line towards
                else
                   disp("initalising new plot")
                   
                    px = rb(i).Position(1);
                    py = rb(i).Position(2);
                    pz = rb(i).Position(3);
                    vx = 0;
                    vy = 0;
                    vz = 0;
                    
                    px2 = posedata.pose.position.x;
                    py2 = posedata.pose.position.y;
                    pz2 = posedata.pose.position.z;
                    vx2 = 0;
                    vy2 = 0;
                    vz2 = 0;
                   
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

%% Write to excel
filename = 'C:\Users\MSI User\Desktop\testdata.xlsx';
% ("Array",filename,~,sheetname,~,range of cells to paste in 'E1:I5'
disp("FINISH")
disp(x)
writematrix(x,filename,'Sheet',1,'Range','B2');
        
        
        