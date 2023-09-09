% SCRIPT_Position_Control
% Position Control for monocopter over optitrack 

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

timer = 0.0;
%% Display data
prb = obj.RigidBody; % previous location of body, used to plot lines
error_past = zeros(1,3);
noudp = [0.0,0.0,0.0]; %no udp matrix
while (ishandle(H))
    if timer > 10.0 && timer < 15.0
        desired = [0.7 0 1.5];
    elseif timer > 15.0 && timer < 20.0
        desired = [0.7 -0.7 1.5];
    elseif timer > 20.0 && timer < 25.0
        desired = [0.0 -0.7 1.5];
    elseif timer > 25.0 && timer < 30.0
        desired = [-0.7 -0.7 1.5];
    elseif timer > 30.0 && timer < 35.0
        desired = [-0.7 0.0 1.5];
    elseif timer > 35.0 && timer < 40.0
        desired = [-0.7 0.7 1.5];
    elseif timer > 40.0 && timer < 45.0
        desired = [0.0 0.7 1.5];
    elseif timer > 45.0 && timer < 50.0
        desired = [0.7 0.7 1.5]; 
    else
        desired = [0.0 0.0 1.5];
    end
%     desired = [0.0 0.0 1.5];
    % Get current rigid body information, this has to be recalled every time for new frame index
    rb = obj.RigidBody;
    
    
    % Output frame information
    %fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);
    %write(up,[0.0,0.0,0.0],"double",computerip,port); % Port oject, data, type of data, computer ip to write to, port
    write(up,noudp,"double",computerip,port); % Port oject, data, type of data, computer ip to write to, port
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
%                 desired = [0.0 0.0 1.0];
                error = desired - feedback;
                derivative_error = error - error_past;
                
                % collective-z
                kp_z = 0.2; %0.6; %0.4; %0.31
                kd_z = 0.8; % 0.6;%0.34; %0.29;  
                input_z = (kp_z * error(3) * 1) + (kd_z * derivative_error(3));
                
                %%% CW GAINS: KP = 0.2, KD = 0.8
                %%% ACW GAINS: KP = 0.2, KD = 1.4

                % cyclic-x
                if eul_deg(2) < 100 %%Gain Scheduled if body roll CW<7.5deg ACW<1.5deg
                    kp_x = 0.25; 
                    kd_x = 0.0625;
                    
                else 
                    %kp_x = 0.15; 
                    %kd_x = 0.45;
                    kp_x = 0.5; 
                    kd_x = 0.125;
                end
%                 input_x = cos(eul_rad(3)) * ((kp_x * error(1) * 1) + (kd_x * derivative_error(1)));
                input_x = ((kp_x * error(1) * 1) + (kd_x * derivative_error(1)));
                
                %%% CW GAINS: KP = 0.5, KD = 0.125
                %%% ACW GAINS : KP = 0.15, KD = 0.0375
                
                % cyclic-y 
                if eul_deg(2) < 100 %%Gain Scheduled if body roll CW <7.5deg ACW<1.5deg
                    kp_y = 0.25; 
                    kd_y = 0.0625;
                else 
                    %kp_y = 0.15; 
                    %kd_y = 0.45;
                    kp_y = 0.5; 
                    kd_y = 0.125;
                end
%                 input_y = sin(eul_rad(3)) * ((kp_y * error(2) * 1) + (kd_y * derivative_error(2)));
                input_y = ((kp_y * error(2) * 1) + (kd_y * derivative_error(2)));

                %%% CW GAINS: KP = 0.5, KD = 0.125
                %%% ACW GAINS : KP = 0.15, KD = 0.0375
                
                input = [input_x, input_y ,input_z];
                
                % fprintf('\t   Position [%f,%f,%f]\n',round(rb(i).Position/1000,2));
                %fprintf('\t   Error [%f,%f,%f]\n', round(error,2));
%                 fprintf('\t   Input [%f,%f,%f]\n', round(input,2));
                fprintf('\t   Timer [%f,%f,%f]\n', round(timer,2));
                fprintf('\t   Desired [%f,%f,%f]\n', round(desired,2));
                fprintf('\t   Bod Att [%f,%f,%f]\n', round(eul_deg,2));
                write(up,input,"double", computerip,port);
                error_past = error;
                noudp = input; %set no udp connection to past input to smooth
                %fprintf('\t   Past error [%f,%f,%f]\n', round(error_past,2));
                x=[x; [now,rb(i).TimeStamp,input_x,input_y,input_z,eul_deg(1),eul_deg(2),eul_deg(3)]];
                timer = timer + (1/30);
            end
        end    
    end
end

%% Write to excel
filename = 'C:\Users\AirLabOptitrack\Desktop\testdata.xlsx';
% ("Array",filename,~,sheetname,~,range of cells to paste in 'E1:I5'
disp("FINISH")
disp(x)
writematrix(x,filename,'Sheet',1,'Range','B2');