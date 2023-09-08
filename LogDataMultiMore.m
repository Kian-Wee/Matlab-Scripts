% SCRIPT_StreamData
%   Stream live data received by from OptiTrack software(Motive)
%   Saves data into excel
%   Some code adapted from https://github.com/kutzer/OptiTrackToolbox

%% Clear workspace, close all figures, clear command window
clear all
close all
clc

%% Variables
import body_obj.*

rate=1/360; %in 1/Hz, how fast the graph updates
bodyname=["moose"]; % multiple bodies allowed
x=["Mtime","Otime","name","x","y","z","qx","qy","qz","qw"]; % Array to store to excel

%% Create OptiTrack object
obj = OptiTrack;
Initialize(obj,'192.168.1.5','multicast'); %Ensure broadcast frame id is on, loop interface is set to this ip and transmission type is set to multicast

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

%% Break loop if keypress to save to excel
DlgH = figure;
H = uicontrol('Style', 'PushButton', ...
                    'String', 'Close', ...
                    'Callback', 'delete(gcbf)');
drawnow

%% Display data
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
        for j = 1:numel(bodyname)
            if strcmp(rb(i).Name,bodyname(j)) == 1

                % If body loses tracking, stop plotting
                if isempty(rb(i).Position)== 0
%                     fprintf('\t %s \n',string(rb(i).Name));
%                     fprintf('\t Position [%f,%f,%f]\n', [round(rb(i).Position(1)/1000,2),round(rb(i).Position(2)/1000,2),round(rb(i).Position(3)/1000,2)]);
%                     fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb(i).Quaternion);
                    my_field = strcat(convertCharsToStrings(rb(i).Name));
                    variable.(my_field).updatepose(rb(i));
%                     disp(variable.(my_field).velocity)
                    x=[x; [now rb(i).TimeStamp string(rb(i).Name) rb(i).Position(1) rb(i).Position(2) rb(i).Position(3) rb(i).Quaternion(1) rb(i).Quaternion(2) rb(i).Quaternion(3) rb(i).Quaternion(4)]];
                end
            end
        end
    end
    drawnow
end

%% Write to excel
filename = 'C:\Users\area_\Documents\testdata.xlsx';
disp("FINISH")
writematrix(x,filename,'Sheet',1,'Range','B2'); % ("Array",filename,~,sheetname,~,range of cells to paste in 'E1:I5')