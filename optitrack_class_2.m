%% Optitrack object
% Modes: "px4" for uXRCE-DDS PX4ROS2 bridge
%        "ros2" to publish as ros2 topic
%        "internal" publish various outputs for use in a matlab script
classdef optitrack_class_2 < handle
   properties
      bodyname="";
      position=[0,0,0];
      quarternion=[0,0,0,0];
      euler=[0,0,0];
      velocity=[0 0 0];
      past_position=[0,0,0];
      past_time=0;
      body_roll;
      roll_rate;
      body_rotation;
   end
   methods
       function init(obj,name)
        obj.bodyname=name;
      end

      function updatepose(obj,rb)
        obj.quarternion = [rb.Quaternion(1) rb.Quaternion(2) rb.Quaternion(3) rb.Quaternion(4)]; % w, x, y, z
        obj.position = [rb.Position(1) rb.Position(2) rb.Position(3)];
        eul = quat2eul(obj.quarternion);
        obj.body_roll=eul
        class(eul)
        isa(eul,'cell')
        a = zeros(3, 1, 'double');
        a = [rad2deg(eul(0)) rad2deg(eul(1)) rad2deg(eul(2))];
        obj.euler = a;
      end

   end
end
