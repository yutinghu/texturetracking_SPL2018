# texturetracking_SPL2018
Feature-based Tracking Method for High-speed Vision System (Version 1.0)

Copyright: 2018/09/24, Yuting Hu
Affiliation: Omni Lab for Intelligent Visual Engineering and Science (OLIVES), 
School of Electrical and Computer Engineering (ECE), Georgia Institute of Technology
Author: Yuting Hu <huyuting@gatech.edu>
Co-authors: Zhiling Long <zhiling.long@ece.gatech.edu>, AlRegib Ghassan <alregib@gatech.edu>

Citation:
If you use the code and data please cite the following in your work:

Y. Hu, Z. Long, and G. AlRegib, "A High-Speed, Real-Time Vision System for Texture Tracking and Thread Counting," 
IEEE Signal Processing Letters, vol. 25, no. 6, pp. 758-762, 2018.  

@article{hu2018high,
  title={A High-Speed, Real-Time Vision System for Texture Tracking and Thread Counting},
  author={Hu, Yuting and Long, Zhiling and AlRegib, Ghassan},
  journal={IEEE Signal Processing Letters},
  volume={25},
  number={6},
  pages={758--762},
  year={2018},
  publisher={IEEE}
}  

1. System requirements

   a. Linux Ubuntu 18.04.1
   b. Processor: Intel(R) Core(TM) i7-4790K CPU @4.00GHz 4.00GHz
   c. Installed memory(RAM):32.0GB
   d. System type: 64-bit Operating System
   e. Matlab version (important): Matlab 2018a  

2. Assumptions

    a. The sampling rate of the high speed camera is high enough to avoid aliasing between two consecutive frames;
    b. The rotation angle between two tested frames is less than 3 degrees;
    c. Two tested frames have enough overlapping regions so that enough feaure point pairs are obtained;
    d. The rotation center for denim fabric is outside the field of view (FOV) of the camera, but close to the FOV;
    e. We have a prior knowledge of the mapping between weave pattern (thread) and the fabric local lattice shown in captrued frames;
    f. The resolution is high enough for blob detection and segmentation;
    g. The lighting is bright and uniform enough for blob detection and segmentation.

3. Main steps of algorithm

    a. Convert RGB image to gray scale
    b. Feature point detection: detectMSERFeatures.m
    c. Feature extraction: extractFeatures_BRISK.m
    d. Feature point matching: matchFeatures.m
    e. Affine matrix estimation: affinematrixestimation.m
    f. Fabric lattice detection: blobdetection.m, globalrepeatdirection.m, and blobtopology.m
    g. Mapping between pixels and local fabric lattices: pixellatticemapping.m
    
4. Expected inputs/outputs

   a. Inputs:
      1) Two test RGB frames (resolution: 256 times 256 times 3) obtained by legacy HSVS with strong lighting: a reference frame and a current frame 
   b. Outputs:
      1) Translation parameters based on fabric basis vectors: A0_coor and B0_coor
      2) Average of detected fabric basis vectors: locallattice_A0 and locallattice_B0
      3) Rotation parameter in cartesian degree: rotationdegree
      4) Translation parameters in pixels: translationinpixels
      5) Computation time in seconds: elapsedtime

5. How to run the code
   
    Example 1:  Run main.m after you finish step 1 and step 2;

    Objective: Calcualte translation parameter and rotation parameter between frame 0 and frame 1
    Ground truth: Rotation degree between frame 0 and frame 1 is 1/6=0.167 degree obtained by munual rotation stage;
    Operation: a. Clear all figures and variables in the memory;
               b. Select two frames you want to compare by specifying image directory and frame number;
               c. Call function featurebasedtracking to compute outputs;
                  for exmample: rotationdegree=0.183 degree

 
    Example 2:  Run main_test.m after you finish step 1 and step 2;
    Objective: Calcualte translation parameter and rotation parameter between frame 1~60 and frame 0
    Ground truth: rotation degree between two consecutive frames is 1/6=0.167 degree;
    Operation: a. Clear all figures and variables in the memory;
               b. Specify totalframe=60;
               c. Call function featurebasedtracking to compute outputs;
