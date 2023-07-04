# Matlab-Scripts
Matlab Scripts for 3D tracking and visualisation with Optitrack and ROS

Use with [ROS Toolbox](https://www.mathworks.com/products/ros.html) and [Optitrack Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/55675-kutzer-optitracktoolbox#:~:text=This%20toolbox%20interfaces%20the%20OptiTrack,(or%20similar)%20software%20package.)



1) Plot live optitrack data in a 3D space

1.1) Stream optitrack data and broadcast over a UDP port

1.2) Stream optitrack data and broadcast to ROS2 topic(See px4_msg installation notes at the bottom

2) Plot live optitrack data and rosdata in a 3D space

2.1) Same as 2 but with velocity and other values

3) Plot Rosbag data

4) Graphs, generate graphs for paper



## Optitrack Tips

### Setup

Often, it is easier and faster to get the coordinate axis right within motive rather than deal with transforms in post. As such, use the calibration square and rotate until it faces the desired orientation . As a plus side, this makes it straightforward for visualisation when running other autonomy scripts too.

After that is done, ensure your body is orientated in the same axis as the defined axis, ie they both face the same north. 


### Marker Setup

The optitrack cameras have a resolution, which means that too small objects cant be tracked when it is far away from the camera, ie too low to the ground when the camera is mounted on a high ceiling. As such, large markers are always prefered to give it the maximum chance of tracking. The downside however is that due to the aforementioned camera resolution, putting large markers close to each other will tend to merge into one incoherent grid of pixels. Also, do note that while motive recommends [4 markers](https://v22.wiki.optitrack.com/index.php?title=Rigid_Body_Tracking), more markers are not always better as the pixels might also be merged. Sometimes less and smaller markers are better, the easiest way to tell this is by toggling between tracking and live camera views and check if there are any pixels between the markers. Dont be afraid to fall below 4 markers if the pixels get merge or the markers gets occluded.

Also, retroreflective sticker markers dont work well on any moving objects as they get obstructed very easily but it might be fine for static landmarks, stick to ball markers.



### Tuning Camera Settings

The objective here is to maximise threshold while minimising exposure. You want to ensure all the balls are detected with a low as an exposure and as high of a threshold as possible to minimise detection of other artifacts as if they might be confused to be a marker as part of a defined object. Though you might feel the urge to tune each camera's exposure and threshold individually, I found that applying same settings to all (by selecting all the camera) was sufficient espically during day/night cycles. As a side note, the camera's LED's can almost always be left on.

Aside from that, removing or blocking any other physical reflective surfaces helps when permissable, else use masks as a last resort espically if the masks are covering the action space as tracking will be lost in those dead spots.

### Network Setup
For the most part, this should be rather straightforward. You can determine the ip address of the optitrack PC by typing ```ipconfig``` into cmd. Ensure that the optitrack is streaming in the data streaming tab in motive and that the ip address is set to the local network and not loopback, and set the z axis to up for convention. If the optitrack stream is still unable to be detected, check that the optitrack pc may be pinged and following that ensure that all the windows firewall rules allow for motive's ports to passthrough.


## Notes on installing custom msgs(PX4_msgs)
Ensure that the submodule are initalised after cloning with git submodule init ; git submodule update after cloning (or use github desktop)

make sure you ros2genmsg() file path refers to root folder ("custom/") that contains px4_msgs ie not inside px4_msgs("custom/px4_msgs") itself.

there is alot of weird errors trying to generate ros2 msgs, need 2 download visual studio to set up the c++ compiler
ros2genmsg(pwd)
errors
Identifying message files in folder 'C:/Users/area_/OneDrive/Documents/GitHub/Matlab-Scripts'..Validating message files in folder 'C:/Users/area_/OneDrive/Documents/GitHub/Matlab-Scripts'..Done.
Done.
[0/1] Generating MATLAB interfaces for custom message packages... 0%Error using ros.internal.ROSProjectBuilder
Current compiler MinGW64 Compiler (C++) is not supported for ROS build. To choose a compiler, run 'mex -setup cpp'.

fix
download visual studio 2022 (this will auto download c++ compilers too i sthink)
mex -setup:'C:\Program Files\MATLAB\R2023a\bin\win64\mexopts\msvcpp2022.xml' C++
mex -setup cpp
