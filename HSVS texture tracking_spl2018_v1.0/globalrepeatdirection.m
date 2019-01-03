%% 5. Lattice detection: blobdetection.m, globalrepeatdirection.m, and blobtopology.m

%%   (2) Use information in the frequency domain to calculate global repeating directions
%%       of repeat pattern: globalrepeatdirection.m

function [globallattice_angle1,globallattice_angle2]...
    =globalrepeatdirection(referenceframe)

%%     Inputs:
%%            referenceframe: reference frame.

%%     Outputs: 
%%            globallattice_direction1 and globallattice_direction2: global
%%            repeat directions of fabric weave pattern in reference frame.

    [rows,~]=size(referenceframe);
    % Calculate 2D-fft of image referenceframe 
    F0=fft2(referenceframe);
    F0_mag=abs(F0);
    F0_mag_enhanced=log(fftshift(F0_mag));
    % We only need one half of FFT domain to get global repeat directions
    F0_mag_enhanced=F0_mag_enhanced(1:rows/2+1,1:end);
    
    % Find peaks in frequency domain, and there are distance limitations on these peaks
    max_loc=zeros(3,2);
    order=1;
    while order<5
        [max_loc(order,1),max_loc(order,2)] = find(F0_mag_enhanced==max(max(F0_mag_enhanced)),1);
        F0_mag_enhanced(max_loc(order,1),max_loc(order,2))=1e-8;
        dist=10000;
        if order>=2;
            for order_tmp=1:(order-1)
           dist_tmp=(max_loc(order,1)-max_loc(order_tmp,1)).^2+(max_loc(order,2)-max_loc(order_tmp,2)).^2;
               if dist_tmp<dist
                   dist=dist_tmp;
               end
            end
            if dist<=9
                max_loc(order,:)=[];
                order=order-1;
            end
        end
        order=order+1;
    end

    % Peaks in shifted frequency domain correspond to the pattern repeat direction in spatial domain
    % There is a minus here since the coorridate is mirror-symmetric with
    % that in image coordinate (blob topology)
    globallattice_angle1=-atand((max_loc(2,2)-max_loc(1,2))/(max_loc(2,1)-max_loc(1,1)));
    globallattice_angle2_tmp=-atand((max_loc(3,2)-max_loc(1,2))/(max_loc(3,1)-max_loc(1,1)));
    globallattice_angle3_tmp=-atand((max_loc(4,2)-max_loc(1,2))/(max_loc(4,1)-max_loc(1,1)));
    
    % Assuming we have a prior knowledge of the angle between two basis
    % vectors of local lattice
    if abs(globallattice_angle2_tmp-globallattice_angle1)>abs(globallattice_angle3_tmp-globallattice_angle1)
        globallattice_angle2=globallattice_angle2_tmp;
    else
        globallattice_angle2=globallattice_angle3_tmp;
    end
    