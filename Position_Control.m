% SCRIPT_StreamData
%   Stream live data received by from OptiTrack software(Motive) to a UDP Link
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
error_past = zeros(1,3);
while (ishandle(H))
    
    % Get current rigid body information, this has to be recalled every time for new frame index
    rb = obj.RigidBody;
    
    
    % Output frame information
    %fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);

    write(up,[0.7,0.7,0.7],"double",computerip,port); % Port oject, data, type of data, computer ip to write to, port
    % Update each rigid body
    for i = 1:numel(rb)

        % Check for correct body
        % This crashes if there is a glitching body, remove any unused body
        % from the asset tab, if not uncomment this line
        if strcmp(rb(i).Name,bodyname) == 1
%         if convertCharsToStrings(rb(i).Name) == convertCharsToStrings(bodyname)
            

                
            % If body loses tracking, stop plotting
            if isempty(rb(i).Position)== 0
                %fprintf(rb(i).Name);
                %fprintf('\t   Position [%f,%f,%f]\n',rb(i).Position/1000);
%                 fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                eul_deg=rad2deg(quat2eul(rb(i).Quaternion,"XYZ"));
                eul_rad= quat2eul(rb(i).Quaternion,"XYZ");
%                 fprintf('\t Degree[x,y,z] [%f,%f,%f]\n',eul_deg);
%                 write(up,[1,round(rb(i).Position(2)/1000,2),round(rb(i).Position(3)/1000,2)],"double",computerip,port); % Port oject, data, type of data, computer ip to write to, port
%                     write(up,[1,0,0],"double",computerip,port); % Port oject, data, type of data, computer ip to write to, port
                
                x_fb = (rb(i).Position(1))/1000;
                y_fb = (rb(i).Position(2))/1000;
                z_fb = (rb(i).Position(3))/1000;
                feedback = [x_fb y_fb z_fb];
                desired = [0.0 0.0 0.5];
                error = desired - feedback;
                derivative_error = error - error_past;
                
                % collective-z
                kp_z = 0.4; %0.31
                kd_z = 0.29;  
                input_z = (kp_z * error(3) * 1) + (kd_z * derivative_error(3));
                
                
                % cyclic-x
                % if eul_deg(2) > 2.6
                kp_x = 0.2; %0.02 , 1.2
                kd_x = 0.015; %0.015 ,4.5
%                 input_x = cos(eul_rad(3)) * ((kp_x * error(1) * 1) + (kd_x * derivative_error(1)));
                input_x = ((kp_x * error(1) * 1) + (kd_x * derivative_error(1)));
                
                
                % cyclic-y 
                kp_y = 0.2;
                kd_y = 0.015;
%                 input_y = sin(eul_rad(3)) * ((kp_y * error(2) * 1) + (kd_y * derivative_error(2)));
                input_y = ((kp_y * error(2) * 1) + (kd_y * derivative_error(2)));
                
                input = [input_x, input_y ,input_z];
                
                % fprintf('\t   Position [%f,%f,%f]\n',round(rb(i).Position/1000,2));
                %fprintf('\t   Error [%f,%f,%f]\n', round(error,2));
                fprintf('\t   Input [%f,%f,%f]\n', round(input,2));
                write(up,input,"double", computerip,port);
                error_past = error;
                %fprintf('\t   Past error [%f,%f,%f]\n', round(error_past,2));
%                 x=[x; [now,rb(i).TimeStamp,input_x,input_y,input_z,eul_deg(1),eul_deg(2),eul_deg(3)]];
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