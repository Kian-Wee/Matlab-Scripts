%% Clear workspace, close all figures, clear command window
clear all
close all
clc

%% Variables
import body_obj.*

rate=1/16; %in 1/Hz, how fast the graph updates
bodyname=["gp"]; % multiple bodies allowed
data_ar=["Mtime","Otime","name","x","y","z","qx","qy","qz","qw","ëuy","ëup","ëur","vx","vy", "vz","pitch_norm"]; % Array to store to excel

%% Create OptiTrack object
obj = OptiTrack;
Initialize(obj,'192.168.1.5','multicast'); %Ensure broadcast frame id is on, loop interface is set to this ip and transmission type is set to multicast

up = udpport("IPV4");
% up = udpport("datagram","OutputDatagramSize",3);
up.EnableBroadcast = true;

% Initalise all bodies at first
rb = obj.RigidBody;
% Create each motive defined rigid body
arr = blanks(numel(rb));

for i = 1:numel(rb)
    my_field = strcat(convertCharsToStrings(rb(i).Name));
    variable.(my_field) = body_obj;
    variable.(my_field).init(convertCharsToStrings(rb(i).Name));
end

disp('Outputing poses for: ');
disp(variable)

computerip="192.168.1.139"; % ip of computer to be written to
port=1234; % port of the computer to be written to
%% Break loop if keypress to save to excel
DlgH = figure;
H = uicontrol('Style', 'PushButton', ...
                    'String', 'Close', ...
                    'Callback', 'delete(gcbf)');
drawnow


exp = ExpAuxiliaryFunctions;
derivatives = exp.circle_setpoints_anti_cw(1,2,2,1); % circle anti_cw setpoints
% derivatives = exp.circle_setpoints_cw(1,-2,2,1); % circle cw setpoints

% vel = load("invert_vel.mat");
% acc = load("invert_acc.mat");
% jer = load("invert_jer.mat");
% sna = load("invert_sna.mat");


%%% step 1 - attitude control:
%%  assign the heading at the period of the time interval which is pi/50 
%   take in measurements, assume the following for now as pseudo
desired_alt = 2;
mea_pos = zeros(3,1); % extract position measurements in real time from opti track 
mea_vel = zeros(3,1); % extract velocity measurements in real time from opti track
mea_acc = zeros(3,1); % extract acceleration measurements in real time from opti track 

mea_pitch = zeros(1,1); % euler angle for disk pitch which is body roll, anyway this is measured at the centre, so no diff
mea_precession_angle = zeros(1,1); % euler angle for disk roll, precession angle
mea_pitch_rate = zeros(1,1); % euler angle for disk pitch rate which is body roll rate
mea_euler = [0,mea_pitch,0]; % default seq is about ZYX


mea_rotation = zeros(1,1); % body yaw angle for azimuth, must be in RAD
mea_xy_pos_mag = zeros(1,1);
mea_xy_vel_mag = zeros(1,1);
trigger = zeros(1,1); % temporary trigger for now to go into offboard mode

% gains
kpos = 0.1;
kvel = 0.1;
kpos_z = 0.1;
prp = [0.1,0.1]; % bodyrate gain
ppq = 0.1; % body acc gain

% init a_des
a_des = zeros(3,1);

% gravity
g = -9.81;

% linear drag coeff
Dx = 0.03;
Dy = 0.03;
Dz = 0.03;
linear_drag_coeff = [Dx,Dy,Dz];

% unit vectors
ey = [0,1,0]; % 1 x 3
ex = [1,0,0];
ez = [0,0,1];

% init error quaternion
error_quat = zeros(1,4);

% init body rates xy
body_rates = zeros(1,2);

% need to insert update rate, loop at 1/time_per_setpt freq which is currently 16 hz
update_rate = derivatives(7,1);

i = 1; % counter
old = 0;
old_precession_rate_angle = 0;
while ishandle(H)
%%
    % Get current rigid body information, this has to be recalled every time for new frame index
    rb = obj.RigidBody;

    % Output frame information
    %fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);

    % Update each rigid body
    for k = 1:numel(rb)

        % Check for correct body
        % This crashes if there is a glitching body, remove any unused body
        % from the asset tab, if not uncomment this line
        for j = 1:numel(bodyname)
            if strcmp(rb(k).Name,bodyname(j)) == 1

                % If body loses tracking, stop plotting
                if isempty(rb(k).Position)== 0
%                     fprintf('\t %s \n',string(rb(i).Name));
%                      fprintf('\t Position [%f,%f,%f]\n', [round(rb(i).Position(1)/1000,2),round(rb(i).Position(2)/1000,2),round(rb(i).Position(3)/1000,2)]);
%                     fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                    my_field = strcat(convertCharsToStrings(rb(k).Name));
                    variable.(my_field).updatepose(rb(k));
                    disp(rad2deg(variable.(my_field).euler));
%                     disp(rad2deg(variable.("gp").euler));
%                    x=[x; [now rb(i).TimeStamp string(rb(i).Name) rb(i).Position(1) rb(i).Position(2) rb(i).Position(3) rb(i).Quaternion(1) rb(i).Quaternion(2) rb(i).Quaternion(3) rb(i).Quaternion(4) variable.(my_field).euler(1) variable.(my_field).euler(2) variable.(my_field).euler(3) variable.(my_field).velocity(1) variable.(my_field).velocity(2) variable.(my_field).velocity(3) variable.(my_field).pitch_norm]];
                end
            end
        end
    end
    drawnow

    mea_pos = transpose(variable.gp.position); % extract position measurements in real time from opti track 
    mea_vel = transpose(variable.gp.velocity); % extract velocity measurements in real time from opti track
    if variable.gp.euler(3) == 0
        mea_pitch = variable.gp.euler(2);
    end

    mea_xy_pos_mag = sqrt((mea_pos(1,:)).^2 + (mea_pos(2,:)).^2);
    mea_xy_vel_mag = sqrt((mea_vel(1,:)).^2 + (mea_vel(2,:)).^2);
%%
    if i > 200
        i = 1;
    end 
    
    %%%% (Precession Rate Tracking)
%     latest_mag = mea_xy_pos_mag - sqrt(8);
%     latest_precession_rate_angle = acos(latest/old);
% 
%     precession_rate = -1 * ((latest-old_mag)/abs((latest-old_mag))) * (latest_precession_rate_angle/update_rate); % if negative in opti-track means correct direction, I multiply by -1 to make it positive so its in the correct direction
% 
%     old_mag = latest_mag;
%     old_precession_rate_angle = latest_precession_rate_angle;
    
    %%%% (Test)
    % xy
%     a_fb_xy = abs((kpos*(derivatives(1,1) - mea_xy_pos_mag)) + (kvel*(derivatives(2,1) - mea_xy_vel_mag))); % xy_magnitude plane since yaw can be easily taken care of 
%     gain = kpos*(derivatives(1,1) - mea_xy_pos_mag)/abs(kpos*(derivatives(1,1) - mea_xy_pos_mag));
%     a_rd = derivatives(2,1) * linear_drag_coeff(1,1);
%     a_des(1,:) = a_fb_xy + derivatives(3,1) - a_rd; % fits into the x axis of ades

    % z (can be used to test, needs to activate hover flaps mode)
%     a_fb_z = kpos_z*(desired_alt - mea_pos(3,1)); % z
%     a_des(3,:) = a_fb_z + g;
%     zd = a_des / norm(a_des); % 3 x 1 = z_desired, if empty it would be 0 0 1  

     % direction 
%    desired_heading = derivatives(6,1);



    %%%% (Actual)
    % xy 
%     a_fb_xy = abs((kpos*(derivatives(1,i) - mea_xy_pos_mag)) + (kvel*(derivatives(2,i) - mea_xy_vel_mag))); % xy_magnitude plane since yaw can be easily taken care of 
%     gain = kpos*(derivatives(1,i) - mea_xy_pos_mag)/abs(kpos*(derivatives(1,i) - mea_xy_pos_mag));
%     a_rd = derivatives(2,i) * linear_drag_coeff(1,1);
%     a_des(1,:) = a_fb_xy + derivatives(3,i) - a_rd; % fits into the x axis of ades 

    % z (can be used to test, needs to activate hover flaps mode)
    a_fb_z = kpos_z*(desired_alt - mea_pos(3,1)); % z
    a_des(3,:) = a_fb_z + g;
    zd = a_des / norm(a_des); % 3 x 1 = z_desired, if empty it would be 0 0 1  

    % direction 
    desired_heading = derivatives(6,i);
    

    %% Bodyrates (for collective thrust test, this entire section can be disabled)
%     qz = eul2quat(mea_euler); % default seq is q = [w x y z]
%     disk_vector = quatrotate(qz,ez); % vector of 1 x 3
%     angle = acos((dot(disk_vector,transpose(zd))/(norm(disk_vector)*norm(zd)))); % will nvr catch up one
%     n = cross(disk_vector,transpose(zd)) / norm(cross(disk_vector,transpose(zd)));
%     B = quatrotate((quatinv(qz)),n);
% 
%     % if angle = 0, it would just be an identity quat matrix
%     error_quat = [cos(angle/2),B*sin(angle/2)];
% 
%     % can alwways break here to make sure shit is running correctly
%     
%     if error_quat(:,1) < 0
%         body_rates = -2 * prp * error_quat(:,2:3);
%     else
%         body_rates = 2 * prp * error_quat(:,2:3);
%     end

    %% Collective thrust (can be used to test)
    cmd_z = dot(transpose(zd),transpose(a_des)); %% command sent to motor, need to include filter to make sure negative cmds dun go thru
    disp("cmd_z: ");
    disp(cmd_z);  

    %% Diff Flat Feedforward component (actual)
%     ff_n = (derivatives(4,i) + linear_drag_coeff(:,1)*derivatives(3,i));
%     v = [derivatives(4,i),0,0]
%     ff_d = cmd_z + linear_drag_coeff(:,1)*dot(zd,v) + linear_drag_coeff(:,1)*dot(ex,v) - linear_drag_coeff(:,3)*dot(zd,v);
%     if ff_d == 0
%         body_rate_ref = 0;
%     else
%         body_rate_ref = ff_n/ff_d;
%     end

    %% testing
%     body_rate_ref = 0;
%     
%     cmd_bodyrate = ppq * (body_rates(:,3) - mea_pitch_rate + body_rate_ref); % gain for cyclic, multiply this to azimuth sin or cos from quadrant, the other value is the desired heading  
% 
%     quadrant = exp.quadrant_output(desired_heading);
%     input = exp.flap_output(mea_rotation,quadrant,gain,desired_heading,abs(cmd_bodyrate));    
%     final_flap_input = input(:,1);
    
    i = i + 1;
%     trigger = trigger + update_rate; % temporary holding


    input = [0,0,cmd_z]; %heading, flap, motor
    fprintf('\t   Input [%f,%f,%f]\n', input);
    write(up,input,"double", computerip,port);


    pause(rate);
end

%% Write to excel
filename = 'C:\Users\area_\OneDrive\Desktop\testdata.xlsx';
disp("FINISH")
writematrix(x,filename,'Sheet',1,'Range','B2'); % ("Array",filename,~,sheetname,~,range of cells to paste in 'E1:I5')
