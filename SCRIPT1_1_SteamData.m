% SCRIPT_StreamData
%   Stream live data received by from OptiTrack software(Motive) to a UDP Link
%   Some code adapted from https://github.com/kutzer/OptiTrackToolbox

%% Clear workspace, close all figures, clear command window
clear all
close all
clc

%% Variables
rate=1/10; %in 1/Hz, how fast the graph updates, set to 0 got instant update
bodyname='test';
computerip="192.168.1.3"; % ip of computer to be written to
port=1234; % port of the computer to be written to

%% Create OptiTrack object, UDP Port Object and initialize
obj = OptiTrack;
Initialize(obj,'192.168.1.5','multicast'); %Ensure broadcast frame id is on, loop interface is set to this ip and transmission type is set to multicast
init = 0;
up = udpport("IPV4");
% up = udpport("datagram","OutputDatagramSize",3);
up.EnableBroadcast = true;

%% Create Live Graph
figure;
%xlim([0 25]);
%ylim([-1.1 1.1]);


%% Display data

% Plot origin
% plot(0,0,'+k');

prb = obj.RigidBody; % previous location of body, used to plot lines

while true
    
    % Get current rigid body information, this has to be recalled every time for new frame index
    rb = obj.RigidBody;
    
    % Output frame information
    fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);

    % Update each rigid body
    for i = 1:numel(rb)

        % Check for correct body
        % This crashes if there is a glitching body, remove any unused body
        % from the asset tab, if not uncomment this line
        if strcmp(rb(i).Name,bodyname) == 1
            
           fprintf(rb(i).Name);
            fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
            fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                
            % If body loses tracking, stop plotting
            if isempty(rb(i).Position)== 0

                if init ~= 0
                    fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
                    fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                    %plot(rb(i).Position(1)/1000,rb(i).Position(2)/1000,'-or');
                    %plot3([rb(i).Position(1)/1000 prb(i).Position(1)/1000],[rb(i).Position(2)/1000 prb(i).Position(2)/1000],[rb(i).Position(3)/1000 prb(i).Position(3)/1000],'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
                    plot3([rb(i).Position(1)/1000 prb(i).Position(1)/1000],[rb(i).Position(2)/1000 prb(i).Position(2)/1000],[rb(i).Position(3)/1000 prb(i).Position(3)/1000], 'b');

                    prb=rb;
                    
                % On initialisation, there is no previous point to plot a line towards
                else
                    disp("initalising new plot")
                    prb=rb; %plot a point
                    plot3([rb(i).Position(1)/1000 prb(i).Position(1)/1000],[rb(i).Position(2)/1000 prb(i).Position(2)/1000],[rb(i).Position(3)/1000 prb(i).Position(3)/1000], '-r');
                    init=1;
                    
                    write(up,[round(rb(i).Position(1)/1000,2),round(rb(i).Position(2)/1000,2),round(rb(i).Position(3)/1000,2)],"uint8",computername,port); % Port oject, data, type of data, computer ip to write to, port
                    write(up,String([round(rb(i).Position(1)/1000,2),round(rb(i).Position(2)/1000,2),round(rb(i).Position(3)/1000,2)]),"uint8",computername,port); % Port oject, data, type of data, computer ip to write to, port
%                     write(up,[1.23,2.34,3.56],"uint8","192.168.1.213",1234);
                    %4*3
                end
%                 %title('Title')
%                 xlabel('x') 
%                 ylabel('y')
%                 zlabel('z') 
%                 grid on;
%                 hold on;
%                 drawnow;
%                 
%                 pause(rate);

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