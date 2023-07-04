%% Create new object
classdef optitrack_class < handle
   properties
      bodyname="";
      body_output_topic="";
      body_pose_pub;
      body_pose_msg;
      mode;
   end
   methods
      function init(obj,obj_name,optitrack_node,mode)
         obj.bodyname=obj_name;
         obj.body_output_topic = strcat(obj_name,"/optitrack_pose");
         obj.mode=mode;
         if mode=="px4" %Assume publishing using uXRCE-DDS PX4ROS2 bridge
            % obj.body_pose_pub = ros2publisher(optitrack_node, "/fmu/in/vehicle_visual_odometry", "px4_msgs/VehicleOdometry","History","keeplast","depth",0,"Reliability","besteffort","Durability","transientlocal");
            obj.body_pose_pub = ros2publisher(optitrack_node, "/fmu/in/vehicle_visual_odometry", "px4_msgs/VehicleOdometry")
            obj.body_pose_msg = ros2message("px4_msgs/VehicleOdometry");
         else 
            obj.body_pose_pub = ros2publisher(optitrack_node, obj.body_output_topic, "geometry_msgs/PoseStamped");
            obj.body_pose_msg = ros2message("geometry_msgs/PoseStamped");
         end
      end

      function sendpose(obj,x,y,z,ww,wx,wy,wz)

      if obj.mode=="px4"
         obj.body_pose_msg.timestamp = uint64((rostime("now").Sec*10^9+rostime("now").Nsec)/1000);
         % obj.body_pose_msg.timestamp_sample = uint64(rostime("now").Sec*10^9+rostime("now").Nsec); %using this causes it not to fuse for some reason
         obj.body_pose_msg.pose_frame = uint8(2); % 1 for NED earth-fixed frame; 2 for FRD arbitrary frame
         obj.body_pose_msg.position(1)=x;
         obj.body_pose_msg.position(2)=-y;
         obj.body_pose_msg.position(3)=-z;
         xx = [[1 0 0];
               [0 cos(pi) -sin(pi)];
               [0 sin(pi) cos(pi)];];
         q1 = rotm2quat(xx);
         qq = quatmultiply(q1,[ww,wx,wy,wz]);
         obj.body_pose_msg.q =single(qq);
      else
         obj.body_pose_msg.header.stamp.sec = cast(posixtime(datetime), "int32");
         obj.body_pose_msg.header.frame_id ='/odom';
         obj.body_pose_msg.pose.position.x = x;
         obj.body_pose_msg.pose.position.y = y;
         obj.body_pose_msg.pose.position.z = z;
         obj.body_pose_msg.pose.orientation.w = ww;
         obj.body_pose_msg.pose.orientation.x = wx;
         obj.body_pose_msg.pose.orientation.y = wy;
         obj.body_pose_msg.pose.orientation.z = wz;

         % % Janky hack for euler on the same topic;
         % eu=quat2eul(ww,wx,wy,wz);
         % obj.body_pose_msg.pose.orientation.w = 0;
         % obj.body_pose_msg.pose.orientation.x = eu(2);
         % obj.body_pose_msg.pose.orientation.y = eu(1);
         % obj.body_pose_msg.pose.orientation.z = eu(0);
         % disp(quat2eul(ww,wx,wy,wz));
      end
        send(obj.body_pose_pub,obj.body_pose_msg); % publish msg
        
      end
   end
end
