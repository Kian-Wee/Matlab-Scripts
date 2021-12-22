% SCRIPT_DisplayData
%   Display data received by from OptiTrack software to the command window.
%
%   M. Kutzer 10Mar2016, USNA

%% Clear workspace, close all figures, clear command window
clear all
close all
clc

%% Create OptiTrack object and initialize
obj = OptiTrack;
Initialize(obj,'192.168.1.25','multicast');

%% Create Live Graph
figure;
%xlim([0 25]);
%ylim([-1.1 1.1]);
x = 1:0.01:25;
y = sin(x);

%% Subscribe to /uav0/mavros/setpoint_position/local to compare
% rosinit('192.168.1.3');
% setpointsub=rossubscriber('/uav0/mavros/setpoint_position/local');



%% Display data
while true
    % Get current rigid body information
    rb = obj.RigidBody;
    % Output frame information
    fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);
    % Update each rigid body
    for i = 1:numel(rb)
        % If body loses tracking, stop plotting
        if rb(i).Name=='BetaVroom' & isempty(rb(i).Position)== 0
            %plot(rb(i).Position(1)/1000,rb(i).Position(2)/1000,'-or');
            plot3(rb(i).Position(1)/1000,rb(i).Position(2)/1000,rb(i).Position(3)/1000,'-or');
            %plot3(xt1,yt1,zt1,xt2,yt2,zt2)
            %plot(0,0,'-hg');
            %plot(x,y,'-hb');
            drawnow;
            hold on;
            % pause(1);
        
            fprintf('- %s, Tracking Status: %d\n',rb(i).Name,rb(i).isTracked);
            if rb(i).isTracked
                fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
                fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
%             else
%                 fprintf('\t Position []\n');
%                 fprintf('\t Quaternion []\n');
            end
        end
    end
end
        
        