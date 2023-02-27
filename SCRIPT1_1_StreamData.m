% SCRIPT_StreamData
%   Stream live data received by from OptiTrack software(Motive) to a UDP Link
%   Also saves data into excel, comment out if not needed
%   Some code adapted from https://github.com/kutzer/OptiTrackToolbox

%reduce exposure increase threshold

%% Clear workspace, close all figures, clear command window
clear all
close all
clc

%% Variables
rate=1/10; %in 1/Hz, how fast the graph updates, set to 0 got instant update
bodyname='DMMC';
computerip="192.168.1.137"; % ip of computer to be written to
port=1234; % port of the computer to be written to
x=["Mtime","Otime","x","y","z","xd","yd","zd"]; % Array to store to excel

%% Create OptiTrack object, UDP Port Object and initialize
obj = OptiTrack;
Initialize(obj,'192.168.1.5','multicast'); %Ensure broadcast frame id is on, loop interface is set to this ip and transmission type is set to multicast
init = 0;
up = udpport("IPV4");
% up = udpport("datagram","OutputDatagramSize",3);
up.EnableBroadcast = true;

%% Break loop if keypress to save to excel
DlgH = figure;
H = uicontrol('Style', 'PushButton', ...
                    'String', 'Break', ...
                    'Callback', 'delete(gcbf)');

%% Display data
prb = obj.RigidBody; % previous location of body, used to plot lines
while (ishandle(H))
    
    % Get current rigid body information, this has to be recalled every time for new frame index
    rb = obj.RigidBody;
    
    
    % Output frame information
    %fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);

    % Update each rigid body
    for i = 1:numel(rb)

        % Check for correct body
        % This crashes if there is a glitching body, remove any unused body
        % from the asset tab, if not uncomment this line
        if strcmp(rb(i).Name,bodyname) == 1
                
            % If body loses tracking, stop plotting
            if isempty(rb(i).Position)== 0
                %fprintf(rb(i).Name);
                %fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
%                 fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                input = [round(rb(i).Position(1)/1000,2),round(rb(i).Position(2)/1000,2),round(rb(i).Position(3)/1000,2)];
                fprintf('\t   Input [%f,%f,%f]\n', input);
                write(up,input,"double", computerip,port);
%                 x=[x; [now,rb(i).TimeStamp,input_x,input_y,input_z,eul_deg(1),eul_deg(2),eul_deg(3)]];
            end
        end    
    end
end

%% Write to excel
filename = 'C:\Users\MSI User\Desktop\testdata.xlsx';
disp("FINISH")
disp(x)
writematrix(x,filename,'Sheet',1,'Range','B2'); % ("Array",filename,~,sheetname,~,range of cells to paste in 'E1:I5')