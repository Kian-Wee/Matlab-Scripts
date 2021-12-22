% SCRIPT_DisplayData
%   Display data received by from OptiTrack software to the command window.
%
%   M. Kutzer 10Mar2016, USNA

%% Clear workspace, close all figures, clear command window
clear all
close all
clc

%% Variables
rate=10; %in Hz, how fast the graph updates
bodyname='BetaVroom';

px =0;
py=0;
pz=0;

%% Create OptiTrack object and initialize
obj = OptiTrack;
Initialize(obj,'192.168.1.25','multicast');
init = 0;

%% Create Live Graph
figure;
%xlim([0 25]);
%ylim([-1.1 1.1]);


%% Display data
while true
    
    % Get current rigid body information
    rb = obj.RigidBody;
    prb = obj.RigidBody; % previous location of body, used to plot lines

    
    % Output frame information
    fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);
    
    % Plot origin
    plot(0,0,'+k');
    
    % Update each rigid body
    for i = 1:numel(rb)

        
        % Check for correct body

        % This crashes if there is a glitching body, remove any unused body
        % from the asset tab
        if rb(i).Name == bodyname 
            
            fprintf(rb(i).Name);
%             fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
%             fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
%             fprintf('\t   Position [%f,%f,%f]\n',prb(i).Position/1000);
%             fprintf('\t Quaternion [%f,%f,%f,%f]\n',prb(i).Quaternion);
                
            % If body loses tracking, stop plotting
            if isempty(rb(i).Position)== 0
                
                if init ~= 0
                    fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
                    fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                    %plot(rb(i).Position(1)/1000,rb(i).Position(2)/1000,'-or');
                    %plot3([rb(i).Position(1)/1000 prb(i).Position(1)/1000],[rb(i).Position(2)/1000 prb(i).Position(2)/1000],[rb(i).Position(3)/1000 prb(i).Position(3)/1000],'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
                    plot3([rb(i).Position(1)/1000 px/1000],[rb(i).Position(2)/1000 py/1000],[rb(i).Position(3)/1000 pz/1000],'b');
                    %plot3([1 10],[1 10],[1 10]);
                    %plot3([1.409585 1.409585],[-2.250146 -2.250146],[0.020222 0.020222]);
                    %p=plot3(rb(i).Position(1)/1000,rb(i).Position(2)/1000,rb(i).Position(3)/1000,'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
                    %plot3(xt1,yt1,zt1,xt2,yt2,zt2)
                    
                    %prb=rb;
                    px = prb(i).Position(1);
                    py = prb(i).Position(2);
                    pz = prb(i).Position(3);
                else
                   disp("initalising new plot")
                   prb=rb;
                   plot3([rb(i).Position(1)/1000 prb(i).Position(1)/1000],[rb(i).Position(2)/1000 prb(i).Position(2)/1000],[rb(i).Position(3)/1000 prb(i).Position(3)/1000], '-r');
                   init=1;
                end
                %title('Title')
                xlabel('x') 
                ylabel('y')
                zlabel('z') 
                grid on;
                hold on;
                drawnow;
                
                %pause(1/rate);

%                 fprintf('- %s, Tracking Status: %d\n',rb(i).Name,rb(i).isTracked);
%                 if rb(i).isTracked
%                     fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
%                     fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
%                 else
%                     fprintf('\t Position []\n');
%                     fprintf('\t Quaternion []\n');
%                 end
            end
        end
    end
end
        
        