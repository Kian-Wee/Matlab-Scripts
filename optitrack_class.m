%% Create new object
classdef optitrack_class < handle
   properties
      bodyname="";
      body_output_topic="";
      body_pose_pub;
      body_pose_msg;
   end
   methods
      function init(obj,obj_name,optitrack_node)
         obj.bodyname=obj_name;
         obj.body_output_topic = strcat(obj_name,"/optitrack_pose");
         obj.body_pose_pub = ros2publisher(optitrack_node, obj.body_output_topic, "geometry_msgs/PoseStamped");
         obj.body_pose_msg = ros2message("geometry_msgs/PoseStamped");
      end

      function sendpose(obj,x,y,z,ww,wx,wy,wz)
        obj.body_pose_msg.pose.position.x = x;
        obj.body_pose_msg.pose.position.y = y;
        obj.body_pose_msg.pose.position.z = z;
%         obj.body_pose_msg.pose.orientation.w = ww;
%         obj.body_pose_msg.pose.orientation.x = wx;
%         obj.body_pose_msg.pose.orientation.y = wy;
%         obj.body_pose_msg.pose.orientation.z = wz;

        % Janky hack for euler on the same topic;
        eu=quat2eul(ww,wx,wy,wz);
        obj.body_pose_msg.pose.orientation.w = 0;
        obj.body_pose_msg.pose.orientation.x = eu(2);
        obj.body_pose_msg.pose.orientation.y = eu(1);
        obj.body_pose_msg.pose.orientation.z = eu(0);
        disp(quat2eul(ww,wx,wy,wz));

        obj.body_pose_msg.header.stamp.sec = cast(posixtime(datetime), "int32");
        obj.body_pose_msg.header.frame_id ='/odom';
        send(obj.body_pose_pub,obj.body_pose_msg); % publish msg
        
      end
   end
end