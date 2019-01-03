%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HSVS feature-based tracking method (Version 1.0) with intermediate results
%% Copyright 2018, Yuting Hu. All rights reserved.
%% Yuting Hu (huyuting@gatech.edu)
%% Run codes in Matlab 2018a

% Citation:
% If you use the code and data please cite the following in your work:
% 
% Y. Hu, Z. Long, and G. AlRegib, "A High-Speed, Real-Time Vision System for Texture Tracking and Thread Counting," 
% IEEE Signal Processing Letters, vol. 25, no. 6, pp. 758-762, 2018.  
% 
% @article{hu2018high,
%   title={A High-Speed, Real-Time Vision System for Texture Tracking and Thread Counting},
%   author={Hu, Yuting and Long, Zhiling and AlRegib, Ghassan},
%   journal={IEEE Signal Processing Letters},
%   volume={25},
%   number={6},
%   pages={758--762},
%   year={2018},
%   publisher={IEEE}
% } 

%% Assupmtion:
%% 1. The sampling rate of the high speed camera is high enough to avoid
%%    aliasing between two consecutive frames;
%% 2. The rotation angle between two tested frames is less than 3 degrees;
%% 3. Two tested frames have enough overlapping regions so that enough feaure
%%    point pairs can be obtained;
%% 4. The rotation center for denim fabric is outside the field of view (FOV)
%%    of the camera, but close to the FOV;
%% 5. We have a prior knowledge of the mapping between weave pattern (thread) and
%%    the fabric local lattice shown in captured frames;
%% 6. The resolution is high enough for blob detection and segmentation;
%% 7. The lighting is bright and uniform enough for blob detection and segmentation.


%% Code objective:
%% Assuming that the sampling rate of the high speed camera is high enough
%% to avoid aliasing between two consecutive frames and assuming that the
%% coordinate origin is in the upper left corner of the image, we can track
%% the movement of back denim fabric surface with blob regions including:
%%      1. Track the translation parameters in local fabric lattice
%%      2. Track the rotation parameters in cartesian degrees
%%  ps: Rotation first, and then translation; the coordinate origin is in
%%      the upper left corner of the image.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% Read a reference frame and a current frame and transfer RGB images to
% grayscale images
referenceframeno=0;
currentframeno=1;
referenceframe = rgb2gray(imread(strcat('images/img_',num2str(referenceframeno),'.bmp'))); 
currentframe = rgb2gray(imread(strcat('images/img_',num2str(currentframeno),'.bmp')));

% Calculate a translation parameter and a rotation parameter between the reference
% frame and the current frame
[locallattice_A0,locallattice_B0,A0_coor,B0_coor,rotationdegree,translationinpixels,elapsedtime]...
    =featurebasedtracking(referenceframe,currentframe);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%