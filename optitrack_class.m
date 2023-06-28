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
            obj.body_pose_pub = ros2publisher(optitrack_node, "/fmu/in/vehicle_mocap_odometry", "px4_msgs/VehicleOdometry","History","keeplast","depth",0,"Reliability","besteffort","Durability","transientlocal");
            obj.body_pose_msg = ros2message("px4_msgs/VehicleOdometry");
         else 
            obj.body_pose_pub = ros2publisher(optitrack_node, obj.body_output_topic, "geometry_msgs/PoseStamped");
            obj.body_pose_msg = ros2message("geometry_msgs/PoseStamped");
         end
      end

      function sendpose(obj,x,y,z,ww,wx,wy,wz)

      if obj.mode=="px4"
         obj.body_pose_msg.timestamp = uint64(rostime("now").Sec*10^9+rostime("now").Nsec);
         obj.body_pose_msg.timestamp_sample = uint64(rostime("now").Sec*10^9+rostime("now").Nsec);
         obj.body_pose_msg.pose_frame_ned= 1; % NED earth-fixed frame;
         obj.body_pose_msg.position(1)=x;
         obj.body_pose_msg.position(2)=y;
         obj.body_pose_msg.position(3)=z;
         obj.body_pose_msg.q(1)=ww;
         obj.body_pose_msg.q(2)=wx;
         obj.body_pose_msg.q(3)=wy;
         obj.body_pose_msg.q(4)=wz;
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