%Pathways to folders containing data
 paths=[];

% edit this to show the paths to the data 

paths{1} =   '\\VAULT2\Vault2Data\Kelson\xi'  
paths{2} =   '\\VAULT2\Vault2Data\Kelson\lambda'   

%File names for various experiments. Psignal files that correspond to commented out expnames will not be processed.

input.expname = [];
 

 input.expname{1} = '18-01-25';
%input.expname{2} = '12-27-17';
%input.expname{3} = '12-27-17';

%Psignal files
psignalfiles=[];

  

input.psignalfiles{1,1} = 'ART_xi_2018_01_25_Phys_2.mat';
input.psignalfiles{2,1} = 'AHL_lambda_2018_01_25_Phys_1.mat';
%input.psignalfiles{1,3} = 'AHL_beta_2017_10_31_Phys_2.mat';
%input.psignalfiles{1,4} = 'AHL_beta_2017_11_01_Phys_1.mat';
%input.psignalfiles{1,5} = 'AHL_beta_2017_11_02_Phys_1.mat';

