function [new_bvals, new_bvecs] = bval_bvec_correction(g,bvals,bvecs)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Reading gradient nonlinearity file and bvals/bvecs
% The I/O below assumes FSL installed (4.0 or higher) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read gradient nonlinearity file
% addpath([getenv('FSLDIR') '/etc/matlab']) 
% filename = 'grad_dev.nii.gz'; 
% g = read_avw(filename); 

% nii1 = load_nii('grad_dev_2slice.nii.gz');
% g = nii1.img;

% Read bvals and bvecs text files 
% bvecs = load('bvecs'); % should be 3xN 
% bvals = load('bvals'); % should be 1xN

% load bvals.txt
% load bvecs.txt

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual correction starts here 
% The following code corrects bvecs and bvals
% for a given voxel (i,j,k) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create matrices 
% L = reshape(squeeze(g(i,j,k,:)),3,3); 

L = reshape(squeeze(g),3,3); 
I = eye(3); 

% correct bvecs and calculate their norm 
v = (I+L)*bvecs; 
n = sqrt(sum(v.^2));

% Normalise corrected bvecs and correct bvals 
new_bvecs = v./repmat(n,3,1) ; 
new_bvals = n.^2.*bvals; 
