%% Rigid Body object
% Used to store variables for a rigid body
classdef body_obj < handle
   properties
      bodyname="";
      position=[0,0,0];
      quarternion=[0,0,0,0];
      euler=[0,0,0];
      velocity=[0 0 0];
      past_position=[0,0,0];
      past_time=0;
      pitch_norm=0;
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
        obj.past_position = obj.position;
        obj.position = [rb.Position(1) rb.Position(2) rb.Position(3)];
        eul = quat2eul(obj.quarternion);
        obj.euler=rad2deg(eul);
        if obj.position-obj.past_position == 0 %Reject 0 division
            obj.velocity=[0 0 0];
        else
            obj.velocity=(obj.position-obj.past_position)/((rb.TimeStamp-obj.past_time)*1000);
        end
        obj.pitch_norm=((1-rb.Position(2)*cos(obj.euler(3)))/sin(obj.euler(3)));
        obj.past_time=rb.TimeStamp;
      end

   end
end
