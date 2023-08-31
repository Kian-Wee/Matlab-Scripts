%% Optitrack object
% Modes: "px4" for uXRCE-DDS PX4ROS2 bridge
%        "ros2" to publish as ros2 topic
%        "internal" publish various outputs for use in a matlab script
classdef optitrack_class_2 < handle
   properties
      bodyname="";
      body_output_topic="";
      body_pose_pub;
      body_pose_msg;
      mode;
      position=[0,0,0];
      quarternion=[0,0,0,0];
      euler=[0,0,0];
      velocity=[];
      body_roll;
      roll_rate;
      body_rotation;
   end
   methods
      function init(obj)
         optiobj = OptiTrack;
         Initialize(optiobj,'192.168.1.5','multicast'); %Ensure broadcast frame id is on, loop interface is set to this ip and transmission type is set to multicast
      end

      function updatepose(obj)
        % Get current rigid body information, this has to be recalled every time for new frame index
        rb = obj.RigidBody;
       
        % Output frame information
        % fprintf('\nFrame Index: %d\n',rb(1).FrameIndex);
        % Update each rigid body
        for i = 1:numel(rb)
            % This crashes if there is a glitching body, remove any unused body
            % from the asset tab, if not uncomment this line
                if isempty(rb(i).Position)== 0 % Check that position array is properly formed(not missing)
                    rb_position = rb(i).Position/1000;
                    rb_quaternion = rb(i).Quaternion; % w, x, y, z
                    obj.position = [rb_position(1),rb_position(2),rb_position(3)];
                    obj
                       
    %                 fprintf('%s\t   Position [%f,%f,%f]\n', rb(i).Name, rb_position);
    %                 fprintf('\t Quaternion [%f,%f,%f,%f]\n',rb_quaternion );
    %                 hydra1_pose_msg.pose.position.x = rb_position(1);
    %                 hydra1_pose_msg.pose.position.y = rb_position(2);            
    %                 hydra1_pose_msg.pose.position.z = rb_position(3);
    %                 hydra1_pose_msg.pose.orientation.w = rb_quaternion(1,1);
    %                 hydra1_pose_msg.pose.orientation.x = rb_quaternion(1,2);
    %                 hydra1_pose_msg.pose.orientation.y = rb_quaternion(1,3);
    %                 hydra1_pose_msg.pose.orientation.z = rb_quaternion(1,4);
    
                    my_field = strcat(convertCharsToStrings(rb(i).Name));
                    variable.(my_field).sendpose(rb_position(1),rb_position(2),rb_position(3),rb_quaternion(1,1),rb_quaternion(1,2),rb_quaternion(1,3),rb_quaternion(1,4));
                end
        end
      end
   end
end
