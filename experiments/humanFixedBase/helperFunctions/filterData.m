function [viconDataFilt,xsensDataFilt ] = filterData(viconDataRaw,xsensDataRaw,filtOpt) 

% filterCapturedData
% Function to pre filter the captured IMU and VICON data using a
% Utilises a Savitzky-Golay FIR filter 
% 
%
%
% Arguments :
% viconDataRaw - raw vicon data structure from capture
% xsensDataRaw - raw xsens data structure from capture
% filtOpt - structure must contain following  : 
%        imu.order - order for filtering imu data
%        imu.window - window for filtering imu data
%        vicon.order - order for filtering vicon data
%        vicon.window - window for filtering vicon data
%
% Author: Naveen Kuppuswamy (naveen.kuppuswamy@iit.it)
% iCub Facility, Istituto Italiano di Tecnologia, 16 February 2016


%% Filter options
xOrder = filtOpt.imu.order;
xWindow = filtOpt.imu.window;

vOrder = filtOpt.vicon.order;
vWindow = filtOpt.vicon.window;

%% xsens
xsensDataFilt.accln = sgolayfilt(xsensDataRaw.accln,xOrder,xWindow) ;
xsensDataFilt.gyro = sgolayfilt(xsensDataRaw.gyro,xOrder,xWindow);
xsensDataFilt.t = xsensDataRaw.t;

%% vicon
viconDataFilt.markers.ltoe = sgolayfilt(viconDataRaw.markers.ltoe,vOrder,vWindow);
viconDataFilt.markers.lhee = sgolayfilt(viconDataRaw.markers.lhee,vOrder,vWindow);
viconDataFilt.markers.lankle = sgolayfilt(viconDataRaw.markers.lankle,vOrder,vWindow);
viconDataFilt.markers.lhip = sgolayfilt(viconDataRaw.markers.lhip,vOrder,vWindow);
viconDataFilt.markers.lsho = sgolayfilt(viconDataRaw.markers.lsho,vOrder,vWindow);
viconDataFilt.markers.rtoe = sgolayfilt(viconDataRaw.markers.rtoe,vOrder,vWindow);
viconDataFilt.markers.rhee = sgolayfilt(viconDataRaw.markers.rhee,vOrder,vWindow);
viconDataFilt.markers.rankle = sgolayfilt(viconDataRaw.markers.rankle,vOrder,vWindow);
viconDataFilt.markers.rhip = sgolayfilt(viconDataRaw.markers.rhip,vOrder,vWindow);
viconDataFilt.markers.rsho = sgolayfilt(viconDataRaw.markers.rsho,vOrder,vWindow);
viconDataFilt.markers.tors = sgolayfilt(viconDataRaw.markers.tors,vOrder,vWindow);
viconDataFilt.markers.imuA = sgolayfilt(viconDataRaw.markers.imuA,vOrder,vWindow);
viconDataFilt.markers.imuB = sgolayfilt(viconDataRaw.markers.imuB,vOrder,vWindow);
viconDataFilt.markers.imuC = sgolayfilt(viconDataRaw.markers.imuC,vOrder,vWindow);
viconDataFilt.analogsMOM = sgolayfilt(viconDataRaw.analogsMOM,vOrder,vWindow);
viconDataFilt.analogsFOR = sgolayfilt(viconDataRaw.analogsFOR,vOrder,vWindow);