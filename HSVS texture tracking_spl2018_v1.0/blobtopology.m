%% 5. Lattice detection: blobdetection.m, globalrepeatdirection.m, and blobtopology.m

%%   (3) Fabric lattice detection: blobtopology.m

function averagelattice=blobtopology(centroids0,globallattice_direction1,...
    globallattice_direction2,blob_template,blobtemplate_mask,f0,flag)

%%     Inputs:
%%            centroids0: the coordinates of all detected MSER centroids in
%%            the reference frame;
%%            globallattice_direction1 and globallattice_direction2: global
%%            repeat directions of fabric weave pattern in reference frame;
%%            referenceframe: the reference frame;
%%            blobtemplate: estimated blob template;
%%            blobtemplate_mask: binary blob template;
%%            flag: indicates whether an MSER feature point in centroids0 is
%%            an individual blob or not.

%%     Outputs: 
%%            averagelattice: average basis vectors of local lattices;


%%     Main steps:
%%     a. Initialzation;
%%     b. FFT-based template matching;
%%     c. Select candidate neighborhood blob centroids for the purpose of lattice
%%        detection based on global directions of repeat patterns;
%%     d. Align all the detected lattices and tick out outliers.

%%    function nnsearch: Copyright, Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     a. Initialzation
    validblobindex=find(flag);
    validbloblistlength=length(validblobindex);
    blob_template_p=(size(blob_template,1)-1)/2;
    [rows,cols]=size(f0);
    % windowsize_factor can be tuned
    windowsize_factor=4;
    detectedlattice=zeros(validbloblistlength,8);
   
% %% Display intermediate results
%    figure,imshow(f0);hold on;
    

for k=1:validbloblistlength
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%     b. FFT-based template matching
    
    % the number 40 can be tuned.
       lattice=zeros(40,2);

    % Here we use inlierpoints0 instead of matchedpoints0
       centroid_tmp=centroids0(validblobindex(k),:); 
    
%     %% Display
%     plot(centroid_tmp(1),centroid_tmp(2),'r*'); 

    regionrow_start=floor(centroid_tmp(2)-windowsize_factor*blob_template_p);
    reginonrow_end=floor(centroid_tmp(2)+windowsize_factor*blob_template_p);
	regioncol_start=floor(centroid_tmp(1)-windowsize_factor*blob_template_p);
	regioncol_end=floor(centroid_tmp(1)+windowsize_factor*blob_template_p);
    centroid_tmp_region=[windowsize_factor*blob_template_p+1,...
        windowsize_factor*blob_template_p+1];
    if regionrow_start<=1
        regionrow_start=1;
        centroid_tmp_region(2)=reginonrow_end-windowsize_factor*blob_template_p;
    end
    if reginonrow_end>=rows
        reginonrow_end=rows;
        centroid_tmp_region(2)=regionrow_start+windowsize_factor*blob_template_p;
    end
    if regioncol_start<=1
        regioncol_start=1;
        centroid_tmp_region(1)=regioncol_end-windowsize_factor*blob_template_p;
    end
    if regioncol_end>=cols
        regioncol_end=cols;
        centroid_tmp_region(1)=regioncol_start+windowsize_factor*blob_template_p;
    end
    f0_region=f0(regionrow_start:reginonrow_end,regioncol_start:regioncol_end);
    f0_region_double=double(f0_region);
    
    % calculate padding
    bx = size(f0_region_double, 2); 
    by = size(f0_region_double, 1);

    % Calculate a correlation map based on 2D-FFT and inverse 2D-FFT
    c = real(ifft2(fft2(f0_region_double) .* fft2(blob_template, by, bx)));

    latticecandidate_count=0;
    
%     % display best match
%     hFig = figure;
%     hAx  = axes;
        
    % find top peak values in the correlation map
    % these peaks correspond to blob centroid candidates
    % this parameter 40 requires tuning
    for j=1:40
        
        [ypeak, xpeak] = find(c == max(c(:)));   
        
        % set the regions correspond to the max peaks into zeros
        setzerosregionrow_start=ypeak-blob_template_p;
        setzerosregionrow_end=ypeak+blob_template_p;
        setzerosregioncol_start=xpeak-blob_template_p;
        setzerosregioncol_end=xpeak+blob_template_p;
        if setzerosregionrow_start<=1
            setzerosregionrow_start=1;
        end
        if setzerosregionrow_end>=by;
            setzerosregionrow_end=by;
        end
        if setzerosregioncol_start<=1
            setzerosregioncol_start=1;
        end
        if setzerosregioncol_end>=bx;
            setzerosregioncol_end=bx;
        end
        
        % lattice candidates are inside the local region of the current centroid
        if latticecandidate_count==0
            lattice(1,:)=[xpeak-blob_template_p,ypeak-blob_template_p];
            latticecandidate_count=latticecandidate_count+1;
            
%             %% display detected matched point
%             plot(lattice(latticecandidate_count,2),...
%             lattice(latticecandidate_count,1),'b*');
        
        else
            lattice_dist=lattice(1:latticecandidate_count,:)...
                -repmat([xpeak-blob_template_p,...
                    ypeak-blob_template_p],latticecandidate_count,1);
            lattice_dist_norm=sum(lattice_dist.^2,2).^0.5;
            % Tuning parameter: 0.6
            % Remove near repetitive peaks within the same blob region
            if min(lattice_dist_norm)>blob_template_p*0.6*sqrt(2);
                latticecandidate_count=latticecandidate_count+1;
                lattice(latticecandidate_count,:)...
                    =[xpeak-blob_template_p,ypeak-blob_template_p];
                
%                 %% display detected matched point
%                 plot(lattice(latticecandidate_count,2),...
%                 lattice(latticecandidate_count,1),'b*');
            
            end
        end

        c(setzerosregionrow_start:setzerosregionrow_end,...
            setzerosregioncol_start:setzerosregioncol_end)...
            =c(setzerosregionrow_start:setzerosregionrow_end,...
            setzerosregioncol_start:setzerosregioncol_end)...
            .*double(~blobtemplate_mask(setzerosregionrow_start...
            +blob_template_p+1-ypeak:setzerosregionrow_end+blob_template_p+1 ...
            -ypeak,setzerosregioncol_start+blob_template_p+1-xpeak...
            :setzerosregioncol_end+blob_template_p+1-xpeak));
        
%         %% Intermediate results: convolution of each local region of inlier
%         %% point and template; zero setting after each max finding iteration
%         figure,imagesc(c);

    end
    lattice(latticecandidate_count+1:end,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%     c. Select candidate neighborhood blob centroids for the purpose of lattice
%%        detection based on global directions of repeat patterns;
    
    % this parameter nearestneighbour_no needs to be tuned
    nearestneighbour_no=12; 
    spatialdist_tmp=lattice-repmat(centroid_tmp_region,latticecandidate_count,1);
    spatialdist_org=sum(spatialdist_tmp.^2,2).^0.5;
    [~,IDX1_tmp]=sort(spatialdist_org);
    IDX1=IDX1_tmp(1:nearestneighbour_no);
    
    % The array of neiborhood centroid candidates
    cellarray1=lattice(IDX1(:),:);
    cellarray1_new=cellarray1(2:end,:)-repmat(cellarray1(1,:),...
    nearestneighbour_no-1,1);
    degree1=atand(cellarray1_new(:,2)./cellarray1_new(:,1)); 
    
    if abs(abs(degree1(1))-abs(globallattice_direction1))<abs(abs(degree1(1))...
            -abs(globallattice_direction2))
        direction2=globallattice_direction2;
    else
        direction2=globallattice_direction1;
    end
    
    % Define a centroid distance based on spatial distance and directional distance
    % spatialdist: spatial distance
      spatialdist=sum(cellarray1_new.^2,2).^0.5;
    % directionaldist: directional distance
      directionaldist=abs(degree1-repmat(direction2,length(cellarray1)-1,1));
    
    % Normalization: the normalizing fatctor can be adjusted
    centroiddist=directionaldist/180.*spatialdist/blob_template_p;
    [~,IDX2_tmp]=sort(centroiddist);
    IDX2 = IDX2_tmp(1:5);
    detectedlattice(k,:)=[cellarray1_new(1,:) cellarray1_new(2,:) ...
        cellarray1_new(IDX2(1),:) cellarray1_new(IDX2(2),:)];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%     d. Align all the detected lattices and tick out outliers
    for n=1:2:3
        if detectedlattice(k,2*n-1)<0
            tmp=detectedlattice(k,2*n-1:2*n);
            detectedlattice(k,2*n-1:2*n)...
                =detectedlattice(k,2*n+1:2*n+2);
            detectedlattice(k,2*n+1:2*n+2)=tmp;
        end
    end
    
    
%     %% Display local lattices
%     figure,imshow(f0_region);hold on;
%     plot(centroid_tmp_region(1),centroid_tmp_region(2),'r*'); hold on;
%     plot(lattice(:,1),lattice(:,2),'b*');
%     
%     hold on;
%     line([centroid_tmp_region(1,1) centroid_tmp_region(1,1)+detectedlattice(k,1)],...
%         [centroid_tmp_region(1,2) centroid_tmp_region(1,2)+detectedlattice(k,2)],'Color','r');
%     line([centroid_tmp_region(1,1) centroid_tmp_region(1,1)+detectedlattice(k,3)],...
%         [centroid_tmp_region(1,2) centroid_tmp_region(1,2)+detectedlattice(k,4)],'Color','g');
%     line([centroid_tmp_region(1,1) centroid_tmp_region(1,1)+detectedlattice(k,5)],...
%         [centroid_tmp_region(1,2) centroid_tmp_region(1,2)+detectedlattice(k,6)],'Color','b');
%     line([centroid_tmp_region(1,1) centroid_tmp_region(1,1)+detectedlattice(k,7)],...
%         [centroid_tmp_region(1,2) centroid_tmp_region(1,2)+detectedlattice(k,8)],'Color','y');
               
end

detectedlattice(validbloblistlength+1:end,:)=[];
latticeoutlier_count=0;
latticeoutlier_idx=zeros(1,validbloblistlength);
for i=1:validbloblistlength
    if ((abs(detectedlattice(i,1)+detectedlattice(i,3))>ceil(blob_template_p/4)+1)...
        ||abs(detectedlattice(i,2)+detectedlattice(i,4))>ceil(blob_template_p/4+1)...
        ||abs(detectedlattice(i,5)+detectedlattice(i,7))>ceil(blob_template_p/4+1)...
        ||abs(detectedlattice(i,6)+detectedlattice(i,8))>ceil(blob_template_p/4+1)...
        ==1)
        latticeoutlier_count=latticeoutlier_count+1;
        latticeoutlier_idx(latticeoutlier_count)=i;
    end
end
latticeoutlier_idx(latticeoutlier_count+1:end)=[];
detectedlattice(latticeoutlier_idx,:)=[];
averagelattice=mean(detectedlattice,1);


% %% Display the average local lattice
% figure,imshow(f0);hold on;
% 
% for k=1:validbloblistlength
% 
%    centroid_tmp=centroids0(validblobindex(k),:); 
%        
%     plot(centroid_tmp(1),centroid_tmp(2),'r*'); 
%     
%     hold on;
%     line([centroid_tmp(1,1) centroid_tmp(1,1)+averagelattice(1)],...
%         [centroid_tmp(1,2) centroid_tmp(1,2)+averagelattice(2)],'Color','r');
%     line([centroid_tmp(1,1) centroid_tmp(1,1)+averagelattice(3)],...
%         [centroid_tmp(1,2) centroid_tmp(1,2)+averagelattice(4)],'Color','g');
%     line([centroid_tmp(1,1) centroid_tmp(1,1)+averagelattice(5)],...
%         [centroid_tmp(1,2) centroid_tmp(1,2)+averagelattice(6)],'Color','b');
%     line([centroid_tmp(1,1) centroid_tmp(1,1)+averagelattice(7)],...
%         [centroid_tmp(1,2) centroid_tmp(1,2)+averagelattice(8)],'Color','y');
%     hold on;
%     
% end
