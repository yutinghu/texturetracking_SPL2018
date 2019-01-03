% Calculate a translation parameter and a rotation parameter between the reference
% frame and the current frame

function [locallattice_A0,locallattice_B0,A0_coor,B0_coor,rotationdegree,...
          translationinpixels,elapsedtime]...
    =featurebasedtracking(referenceframe,currentframe)

%% Inputs and outputs of function featurebasedtracking:

%% Inputs:
%% (1) referenceframe and currentframe: two test frames in grayscale

%% Outputs:
%% (1) A0_coor and B0_coor: translation parameters projected onto fabric basis vectors;
%% (2) locallattice_A0 and locallattice_B0: detected fabric basis vectors;
%% (3) rotationdegree: rotation parameter in cartesian degree;
%% (4) translationinpixels: translation parameters in pixels;
%% (5) elapsedtime: computation time in seconds.


%% Main steps of the feature-based tracking algorithm:
%% 1. Feature point detection: detectMSERFeatures.m
%%    Maximally Stable Extremal Regions (MSER)
%% 2. Feature extraction: extractFeatures_BRISK.m
%%    Binary Robust Invariant Scalable Keypoints (BRISK)
%% 3. Feature point matching: matchFeatures.m
%% 4. Affine matrix estimation: affinematrixestimation.m
%% 5. Lattice detection: blobdetection.m, globalrepeatdirection.m, and blobtopology.m
%% 6. Mapping between pixels and local fabric lattices: pixellatticemapping.m

%% function detectMSERFeatures and matchFeatures: Copyright, Matlab

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Feature point detection: detectMSERFeatures.m
%%     Detect MSER features in a reference frame and a current frame;

%%     Inputs: 
%%             referenceframe and currentframe: two test frames in grayscale

%%     Outputs: 
%%             regions0 and regions1: MSERRegions objects, containing
%%             information of detected MSER features, such as blob centroids
%%             and pixel lists;
%%             centroids0: the coordniates of all MSER centroids in
%%             referenceframe;
%%             MaxAreaVariation: maximum area variation between extremal
%%             regions; its typical values range from 0.1 to 1.0; increasing
%%             this value returns more regions, but there regions may be less stable.

    regions0 = detectMSERFeatures(referenceframe,'MaxAreaVariation',0.2);
    centroids0 = cat(1, regions0.Location);
    regions1 = detectMSERFeatures(currentframe,'MaxAreaVariation',0.2);
    
% %% Intermediate results for MSER features
%     figure; imshow(referenceframe); hold on;
%     plot(regions0); % by default, plot displays ellipses and centroids
%     figure; imshow(currentframe); hold on;
%     plot(regions1); % by default, plot displays ellipses and centroids

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. Feature extraction: extractFeatures_BRISK.m
%%     Extract a BRISK descriptor for each MSER region;

%%     Inputs:
%%           regions0 and regions1: MSER regions

%%     Outputs:
%%           features0 and features1: BRISK feature vectors/descriptors for
%%           each valid MSER feature points;
%%           valid_points0 and valid_points1: valid MSER feature points corresponding
%%           to feature vectors;"valid" means excluding the feature points in
%%           the edges of images;

%% PS: Revise Matlab built-in function extractFeatures;
%% Add another return varible to track the indices of valid feature points: valid_indices

    [features0, valid_points0] = extractFeatures(referenceframe,...
        regions0,'Method','BRISK');
    [features1, valid_points1] = extractFeatures(currentframe,...
        regions1,'Method','BRISK');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. Feature point matching: matchFeatures.m

%%     Intputs:
%%            features0 and features1: BRISK feature vectors in the refernece frame
%%            and the current frame;

%%     Outputs:
%%            indexPairs: index pairs for matched feature points in valid_points0
%%            and valid_points1;
%%            matchedPoints0 and matchedPoints1: coordinates of matched pairs.
%% PS: Call Matlab built-in function matchFeatures;

    indexPairs = matchFeatures(features0, features1);
    matchedPoints0 = valid_points0(indexPairs(:, 1), :);
    matchedPoints1 = valid_points1(indexPairs(:, 2), :);
    
% %% Intermediate results for matched MSER features
%     figure; 
%     showMatchedFeatures(referenceframe,currentframe,matchedPoints0,...
%         matchedPoints1,'false');
%     title('Candidate point matches');
%     legend('Matched points 0','Matched points 1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4. Affine matrix estimation: affinematrixestimation.m

%%     Inputs:
%%            matchedPoints0 and matchedPoints1: coordinates of matched pairs;

%%     Outputs: 
%%            tform: estimated affine matrix from matched feature point pairs;
%%            inliers: logical inlier indices corresponding to indexPairs; "1"
%%            in inliers means inlier feature point; "0" notes an outlier;
%%            indexPairs_updated: inlier indices in centroid lists (centroids0
%%            and centroids1);
%%            rotationdegree: rotation parameters in cartesian degrees;
%%            translationinpixels: translation parameters in pixels between
%%            the reference frame and the current frame.

%% PS: Revise Matlab built-in function extractFeatures;
%% Add another return varible to track indices for inlier feature points: inliers

    [tform,~]= ...
       affinematrixestimation(matchedPoints0,matchedPoints1,'affine',...
       'Confidence',99.8,'MaxDistance',1);


% %% Display intermediate results
% %% Show inlier MSER feature point pairs
%    figure; showMatchedFeatures(referenceframe,currentframe,...
%        centroids0(indexPairs_updated(:,1),:),...
%        centroids1(indexPairs_updated(:,2),:),'false');
%    title('Inlier point matches');
%    legend('Inlier points 0','Inlier points 1');

   % rotationdegree: a rotation paramter in a cartesian degree
   rotationdegree=(-asind(tform.T(1,2))+asind(tform.T(2,1)))/2; 
   translationinpixels=tform.T(3,1:2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5. Lattice detection: blobdetection.m, globalrepeatdirection.m, and blobtopology.m

%%   (1) Individual blob detection based on detected MSER feature points:
%%       blobdetection.m;

%%     Inputs:
%%            referenceframe: a reference frame;
%%            centroids0: the coordinates of all detected MSER centroids in
%%            the reference frame;
%%            regions0: MSERRegions objects, containing information of 
%%            detected MSER features in the reference image.

%%     Outputs: 
%%            blobtemplate: estimated blob template;
%%            blobtemplate_mask: binary mask of the blob template;
%%            flag: indicates whether an MSER feature point in centroids0 is
%%            an individual blob or not;

      [blobtemplate,blobtemplate_mask,flag]...
          =blobdetection(referenceframe, centroids0,regions0);
    
%%   (2) Use information in the frequency domain to calculate global repeating directions
%%       based on placement of repeat pattern: globalrepeatdirection.m

%%     Inputs:
%%            referenceframe: a reference frame.

%%     Outputs: 
%%            globallattice_direction1,globallattice_direction2:
%%            global repeat directions of fabric weave pattern in the reference frame.

        [globallattice_direction1,globallattice_direction2]...
            =globalrepeatdirection(referenceframe);

%%   (3) Fabric lattice detection: blobtopology.m

%%     Inputs:
%%            centroids0: the coordinates of all detected MSER centroids in
%%            the reference frame;
%%            globallattice_direction1 and globallattice_direction2: global
%%            repeat directions of fabric weave pattern in the reference frame;
%%            referenceframe: the reference frame;
%%            blobtemplate: estimated blob template;
%%            blobtemplate_mask: binary blob template;
%%            flag: indicates whether an MSER feature point in centroids0 is
%%            an individual blob or not.

%%     Outputs: 
%%            averagelattice: average basis vectors of local lattices;

        averagelattice=blobtopology(centroids0,globallattice_direction1,...
            globallattice_direction2,blobtemplate,blobtemplate_mask,...
            referenceframe,flag);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
%% 6. Mapping between pixels and local fabric lattices: pixellatticemapping.m

%%     Inputs:
%%           tform: geometric affine matrix estimated from matched feature point pairs;
%%           averagelattice: basis vectors of average local lattices;


%%     Outputs: 
%%            locallattice_A0 and locallattice_B0: basis vectors of detected local lattice;
%%            A0_coor and B0_coor: translation parameters based on fabric basis vectors;
%%            locallattice_A0 and locallattice_B0: avarage of detected local lattice;
%%            A0_coor and B0_coor: translation parameters based on fabric basis vectors;
%%            elapsedtime: computation time in seconds.

    [locallattice_A0,locallattice_B0,A0_coor,B0_coor]= ...
        pixellatticemapping(tform,averagelattice);
    
    elapsedtime=toc; % average compuation time in seconds
         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


