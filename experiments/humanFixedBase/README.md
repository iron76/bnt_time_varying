Experiments related to the human dynamics estimation project
============================================================

This repo presents the results on Maximum A Posteriori estimation of the dynamics of human subjects performing a simple Bowing task. From the physical characteristics of the subject, a simple 3 link, 2 DoF URDF is built and utilised for the computations. The subjects are equipped with a IMU placed on their chest and their movements and contact forces with the ground are captured using a VICON motion capture system. The data from the VICON is assumed to be preprocessed into the format of a matlab .mat file using appropriate libraries (see https://github.com/claudia-lat/VICON-C3D-Analysis for the compatible C3D based analysis of Vicon data).


Requirements :

1. Matlab 2012 or above (not tested on lower version).
2. You need to have properly installed BNT_time_varying (refer to https://github.com/iron76/bnt_time_varying)

Data characteristics
--------------------

The data required for the analysis must be independently downloaded into a subfolder called data. The data consists of 2 files like these: 

1. VICONsaveData.mat
2. imuExtractedData.mat

Subject information
The subject specific URDF must also be located into a subfolder called data

Execution
---------

The execution sequence is as follows : 

1. synchroniseCaptureData
2. organiseSensorData
3. computeSubjectSpecificURDFParams (createUrdfModelFromSubjectParams: step to be computed only for the first time)
4. loadModelFromURDF
5. computeLinkSensorFrame
6. organiseBERDYCompatibleSensorData
7. main


The output of 1,2,5,6,7, produces a .mat file in the subfolder called intermediateDataFiles)

Tests
-------------------
There are some tests available in the tests subfolder. Add this folder to path to execute the tests. For example, predictSensorMeasurement computes an RNEA based forward dynamics and then predicts sensor measurements based on the computed dynamics.

Frame Information
-----------------
For all of the compuations presented here, the frames are assigned as follows : 

![Frame assignment](https://github.com/iron76/bnt_time_varying/blob/dev/naveen/experiments/humanFixedBase/data/framesViconBowingTaskExperiment.jpg)

![Link Frames assignment](https://github.com/iron76/bnt_time_varying/blob/dev/naveen/experiments/humanFixedBase/data/frames_diagram.png)




