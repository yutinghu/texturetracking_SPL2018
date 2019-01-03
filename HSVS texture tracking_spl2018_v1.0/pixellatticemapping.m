%% 6. Mapping between pixels and local fabric lattices: pixellatticemapping.m

function [locallattice_A0,locallattice_B0,A0_coor,B0_coor]...
    =pixellatticemapping(tform,averagelattice)
    
%%     Inputs:
%%            tform: geometric affine matrix estimated from matched feature point pairs;
%%            averagelattice: basis vectors of average local lattices;

%%     Outputs: 
%%            locallattice_A0 and locallattice_B0: basis vectors of detected local lattice;
%%            A0_coor and B0_coor: translation parameters based on fabric basis vectors;
%%            elapsedtime: computation time in seconds.

%%     Main steps:
%%     a. Initialization
%%     b. Project the translational motion vector onto local lattice
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     a. Initialization

% Assumption of origin: upper left corner (0,0) ; rotation first and then
% translation

% Read tranlation motion vector from tform;
    translationmotionvector_fabric=tform.T(3,1:2);
% Read basis vectors of local lattice;
    unitvectors_cand=[averagelattice(1:2);averagelattice(3:4);...
        averagelattice(5:6);averagelattice(7:8);];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%     b. Project the translational motion vector onto local lattice
            
% Care about positive values in iFabricthread and corresponding unit basis vectors
    iFabricthread=unitvectors_cand*translationmotionvector_fabric'...
        ./norm(unitvectors_cand)./norm(translationmotionvector_fabric);
    iFabricthread_positive=find(iFabricthread>0);  

    if length(iFabricthread_positive)==1
        iFabricthread_dismiss=find(iFabricthread==min(iFabricthread));

        % Assume the motion vector point is within the correct local
        % lattice triangular;
        % Select two local lattices from four based on the
        % direction of motion vector

        for i=1:length(iFabricthread)
            if i~=iFabricthread_dismiss && i~=iFabricthread_positive
               unitvectors_cand_tmp1=unitvectors_cand(iFabricthread_positive,:);
               unitvectors_cand_tmp2=unitvectors_cand(i,:);
               S_total=det([1 0 0; 1 unitvectors_cand_tmp1(1) unitvectors_cand_tmp1(2);...
                   1 unitvectors_cand_tmp2(1) unitvectors_cand_tmp2(2)])/2;
               S_1=det([1 0 0; 1 unitvectors_cand_tmp1(1) unitvectors_cand_tmp1(2);...
                   1 translationmotionvector_fabric(1) translationmotionvector_fabric(2)])/2;
               S_2=det([1 0 0; 1 translationmotionvector_fabric(1) ...
                   translationmotionvector_fabric(2);
                   1 unitvectors_cand_tmp2(1) unitvectors_cand_tmp2(2)])/2;
               S_3=det([1 translationmotionvector_fabric(1) translationmotionvector_fabric(2);...
                   1 unitvectors_cand_tmp1(1) unitvectors_cand_tmp1(2);...
                   1 unitvectors_cand_tmp2(1) unitvectors_cand_tmp2(2)])/2;

                if abs(S_total-S_1-S_2-S_3)<1e-3
                    iFabricthread_positive(2)=i;
                end
            end
            if length(iFabricthread_positive)==2
                break
            end
        end
    end
           
   unitvectors_norm=[norm(unitvectors_cand(iFabricthread_positive(1),:));...
       norm(unitvectors_cand(iFabricthread_positive(2),:))];

   % Basis vectors of local lattice: A0 and B0
   % A0: longer basis vector, B0: shorter one
    if unitvectors_norm(1)<unitvectors_norm(2)
        locallattice_A0=unitvectors_cand(iFabricthread_positive(2),:);
        locallattice_B0=unitvectors_cand(iFabricthread_positive(1),:);

    else
        locallattice_B0=unitvectors_cand(iFabricthread_positive(2),:);
        locallattice_A0=unitvectors_cand(iFabricthread_positive(1),:);
    end

    % Project pixel-based translation parameters onto lattice-based
    % translation parameters
    A0_coor=locallattice_A0*translationmotionvector_fabric'./(unitvectors_norm(2).^2);
    B0_coor=locallattice_B0*translationmotionvector_fabric'./(unitvectors_norm(1).^2);   