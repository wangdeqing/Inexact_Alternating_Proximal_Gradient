% Author: Deqing Wang
% Email: deqing.wang@foxmail.com
% Website: http://deqing.me/
% Affiliation: Dalian University of Technology, China
%              University of Jyväskylä, Finland
% Date: April 22, 2021
% Desctirption: Tensor decomposition of an ERP dataset using the
% inexact alternating proximal gradient (iAPG) algorithm.
%
%
% Citation Information:
% D. Wang and F. Cong, An inexact alternating proximal gradient algorithm
% for nonnegative CP tensor decomposition,
% Science China Technological Sciences, 2021. Accepted.
%
%
%% Please install the following toolboxes before running this example
%
% Tensor Toolbox
% http://www.tensortoolbox.org/
%
% EEGLAB
% https://sccn.ucsd.edu/eeglab/download.php
%
%

%%
clear;
close all

%% Data storage path
ERP_data_path = [pwd filesep 'Tutorialdataset2'];
if ~exist(ERP_data_path,'dir')
    % Data can be downloaded manually from http://www.erpwavelab.org/.
    mkdir(ERP_data_path);
    fprintf('Downloading ERP data (115 MB) ......\n');
    data_url = 'http://www.imm.dtu.dk/~mm/downloads/Tutorialdataset2.zip';
    data_filename = 'Tutorialdataset2.zip';
    data_file=websave(data_filename,data_url);
    unzip(data_file,ERP_data_path);
    fprintf('Data has been downloaded and unzipped to\n\n    ''%s''\n\n',...
        ERP_data_path);
end

%% Tensor Representation of the data
% Explanations of the tensor decomposition of this dataset can be found in
% Wang et al., 2018, https://doi.org/10.1016/j.jneumeth.2018.07.020
%
% Tensor parameters
ChannelMode   = 64;
FrequencyMode = 61;
TimeMode      = 72;
SubjectMode   = 14;
ConditionMode = 2;
% Allocate tensor space for all ERP data.
% Fifth-order tensor
DataTensor=zeros(ChannelMode,FrequencyMode,TimeMode,SubjectMode,ConditionMode);
% Load each ERP data to one tensor
dir_items = dir([ERP_data_path filesep '*.mat']);
for FileIndex=1:length(dir_items)
    % File storaging position
    FileName = dir_items(FileIndex).name;
    FilePath = dir_items(FileIndex).folder;
    CurrentDataPath = fullfile(FilePath,FileName);
    % Load data
    CurrentData = load(CurrentDataPath);
    ConditionSign = FileName(1);
    if strcmp(ConditionSign,'L')
        DataTensor(:,:,:,FileIndex,1) = abs(CurrentData.WT);
    end
    if strcmp(ConditionSign,'R')
        DataTensor(:,:,:,FileIndex-14,2) = abs(CurrentData.WT);
    end
end
% Fourth-order tensor
DataTensor4th=reshape(DataTensor,ChannelMode,FrequencyMode,TimeMode,...
    SubjectMode*ConditionMode);
% Third-order tensor
DataTensor3rd=reshape(DataTensor,ChannelMode,FrequencyMode*TimeMode,...
    SubjectMode*ConditionMode);

%% Preparation of tensor decomposition
TensorTrue = tensor(double(DataTensor3rd)); % The third-order tensor
N = ndims(TensorTrue);

% Tensor Decomposition Parameters
R = 20; % The pre-defined number of components

%
ModeSizes = size(TensorTrue);
FR = zeros(N,1); % FR is the ratio of the number of data entries to the degrees of freedom
inv_FR = zeros(N,1);
for ii = 1:N
    FR(ii,1) = prod(ModeSizes) / (R * (ModeSizes(ii) + prod(ModeSizes([1:ii-1 ii+1:end])) - R));
    inv_FR(ii,1) = 1/FR(ii,1);
end
IndicatorDL = nthroot(prod(inv_FR),N); % DL: difficulty level
sIndicatorDL = sprintf('%.3f',IndicatorDL);
fprintf(['ModeSize:\t' mat2str(ModeSizes) '\n' 'RankSize: \t' mat2str(R) '\n']);
fprintf(['DL:\t\t\t' sIndicatorDL '\n']);

% Empirical rule for selecting the number of inner iterations
if IndicatorDL < 0.1
    J = 10;
elseif IndicatorDL >= 0.1
    J = 20;
end
fprintf(['The number of inner iterations is ' mat2str(J) '.\n']);
InnerIter_v = J*ones(N,1);

%% Start of the NCP tensor decomposition using iAPG algorithm
rng('shuffle');
[A,Out] = ncp_iapg(TensorTrue,R,'maxiters',99999,'tol',1e-6,...
    'init','random','printitn',1,'inner_iter',InnerIter_v,...
    'maxtime',1200,'stop',2,'printitn',1);

%%
fprintf('Elapsed time is %4.6f seconds.\n',Out.time(end));
fprintf('Solution relative error = %4.4f\n\n',Out.relerr(2,end));

%%

