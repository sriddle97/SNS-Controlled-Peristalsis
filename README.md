# SNS_Controlled_Peristalsis
This repository contains code that runs the SNS coupled oscillator simulation for control of the peristaltic locomotion of a worm robot as well as the stability analysis of the system's limit cycle. The 2D purely kinematic model can be found in the "Controller Design and 2D Simulation" folder. The 3D model which incorporates soft bodied physics and dynamics is located in the "3D Soft Physics Model" folder. This model much more accurately simulates the behavior of a physical worm robot.

## 2D Simulation
### Running the model
This can be found in SNS_Worm_Sim.m. Everything you need to run the model is contained in this file.

### Stability Analysis
This is contained in the files labeled stability.m and Riddle_cycle_new.m. Both must be located in the same folder for stability.m to run properly.

## 3D Simulation
All files needed to build a worm robot model and run the simulation are contained in the 3D Soft Physics Model folder.

### Building the Model
Mujoco reads models into the simulation environment from an xml file. The xml file for this robot can be auto-generated using the master_worm.m file in the worm_modeling folder. Please note that changes to the model will necessitate minor changes to the code (accounting for new sensor length tolerances, more segments, etc.).

### Running the Model
The python script 3D_SNS_Worm_jupyter_working.ipynb runs the simulation using the model found in worm_3_seg_working.xml. Please note you may need to change the file path mujoco uses to read the xml file in the code.

### Necessary Packages
To run this model you will need to have the both mujoco and the sns_toolbox package installed. The links below will take you to the necessary pages:

https://github.com/wnourse05/SNS-Toolbox.git

https://github.com/deepmind/mujoco/releases