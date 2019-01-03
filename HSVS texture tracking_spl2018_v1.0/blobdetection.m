%% 5. Lattice detection: blobdetection.m, globalrepeatdirection.m, and blobtopology.m

%%   (1) Individual blob detection based on detected MSER feature points:
%%       blobdetection.m;

function [blobtemplate,blobtemplate_mask,flag]=blobdetection(f0,centroids,regions)

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

%%     Main steps:
%%     a. Initialzation
%%     b. Collect blob areas from detected valid MSER regions and caculate a
%%        histogram of blob areas and the mean value of blob areas;
%%     c. Select individual blobs based on a blob area threshold;
%%     d. Learn a blob template based on detected individual blobs;
%%     e. Tick out single blob outliers.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     a. Initialzation

   % flag: indicates whether an MSER feature point in centroids0 is an individual
   % blob or not;
   flag=zeros(length(centroids),1);
   % calculate the row and column size of the input image
   [f0_row,f0_col]=size(f0);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     b. Collect blob areas from detected valid MSER regions and caculate a 
%%        histogram of blob areas and the mean value of blob areas;
    valid_indices_length=length(regions);
    blobarea=zeros(1,valid_indices_length);
    blobmeanintensity=zeros(1,valid_indices_length);
    for k = 1:valid_indices_length
        currentpixellist=regions(k).PixelList;
        blobarea(k)=length(currentpixellist);
        idx_tmp=sub2ind([f0_row,f0_col], currentpixellist(:,2), currentpixellist(:,1));
        sum_tmps=sum(f0(idx_tmp));
        blobmeanintensity(k)=sum_tmps/blobarea(k);
        clear currentpixellist
    end
    
%     %% Intermediate results
%     figure,hist(blobmeanintensity,50);

    blobmeanintensity_mean=mean(blobmeanintensity);
    % This paramter 250 requires tuning.
    [hist_N,hist_X]=hist(blobarea,250);
    % Assuming that the biggest bin corresponds to single blobs
    [~,hist_N_max_idx] = max(hist_N);
    blobarea_threshold=hist_X(hist_N_max_idx);
    
%     %% Intermediate results for histogram of areas of MSER feature points 
%     figure,hist(blobarea,50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     c. Select individual blob based on a blob area threshold;

    % this parameter is odd and it requires tuning
    blobtemplate_row=floor(sqrt(blobarea_threshold)*2*1.2);   
    blobtemplate_col=blobtemplate_row;
    if mod(blobtemplate_row,2) == 0
        blobtemplate_row=blobtemplate_row+1;
        blobtemplate_col=blobtemplate_col+1;
    end
    blobtemplate_p=(blobtemplate_row-1)/2;
    blobtemplate=zeros(blobtemplate_row,blobtemplate_col);
    
% %%     Display
%     figure; imshow(f0); hold on;
    
    count_tmp=0;
    normalarea_loc=zeros(1,valid_indices_length);
    % find reasonable area in bright regions
    for k=1:valid_indices_length
        % this parameter 0.3 and 0.8 need to be tuned
        if (abs(blobarea(k)-blobarea_threshold)/blobarea_threshold<0.2)...
                &&(blobmeanintensity(k)>blobmeanintensity_mean*0.8)...
                &&(sum(floor(regions(k).Location)...
                    -[blobtemplate_p blobtemplate_p]>=[0 0])==2)...
                &&(sum(floor(regions(k).Location)...
                    +[blobtemplate_p blobtemplate_p]<=[f0_col f0_row])==2)
           count_tmp=count_tmp+1;
           normalarea_loc(count_tmp)=k;
           flag(k)=1;
           
% %%            Display
%            plot(regions(k)); 
%            plot(regions(k),'showPixelList',true,...
%                'showEllipses', false); 
%            hold on;

        end
    end
    normalarea_loc(count_tmp+1:end)=[];
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     d. Learn a blob template based on detected individual blobs;
    
    candidate=zeros(blobtemplate_row,blobtemplate_col,length(normalarea_loc));
    normalarea_loc_length=length(normalarea_loc);
    for k = 1:length(normalarea_loc)
        tmp=normalarea_loc(k);
        centroid_tmp=floor(centroids(tmp,:));
        clear currentpixellist
        mask=zeros(blobtemplate_row,blobtemplate_col);
        currentpixellist=regions(tmp).PixelList;
        % align current mask with pixel lists of blob template candidates
        % in original image
        for j=1:length(currentpixellist)
            if double(currentpixellist(j,1))+blobtemplate_p+1-double(...
                    centroid_tmp(1,1))>0
                mask(double(currentpixellist(j,2))+blobtemplate_p+1 ...
                -double(centroid_tmp(1,2)),double(currentpixellist(j,1))...
                    +blobtemplate_p+1-double(centroid_tmp(1,1)))=1;
            end
        end
       
        if centroid_tmp(2)-blobtemplate_p>0 && centroid_tmp(1)-blobtemplate_p>0 ...
                && size(mask,2)==size(candidate,2)
            candidate(:,:,k)=f0(centroid_tmp(2)-blobtemplate_p:centroid_tmp(2)...
                +blobtemplate_p,centroid_tmp(1)-blobtemplate_p:centroid_tmp(1)...
                +blobtemplate_p);
            candidate(:,:,k)=candidate(:,:,k).*double(mask);

            blobtemplate=blobtemplate+candidate(:,:,k);
        else
            normalarea_loc_length=normalarea_loc_length-1;
        end
    end
    clear tmp
    % blobtemplate: learned blob template by averaging individual blob candidates
    blobtemplate=blobtemplate/length(normalarea_loc);
    blobtemplate_mask=bradley(blobtemplate);
    blobtemplate=blobtemplate.*blobtemplate_mask;
    
% %% Intermediate results: learned blob template
%     figure,imshow(blobtemplate,[]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     e. Tick out single blob outliers.
    index=find(flag==1);
    orientations = cat(1, regions.Orientation);
    orientations_tmp=orientations(index);
    index_tmp=(orientations_tmp>mean(orientations_tmp)+std(orientations_tmp));
    flag(index(index_tmp))=0;
