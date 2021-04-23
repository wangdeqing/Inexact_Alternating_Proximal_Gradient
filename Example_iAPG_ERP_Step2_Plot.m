% Author: Deqing Wang
% Email: deqing.wang@foxmail.com
% Website: http://deqing.me/
% Affiliation: Dalian University of Technology, China
%              University of Jyväskylä, Finland
% Date: April 22, 2021
% Desctirption: Plot ERP components extracted by NCP using iAPG algorithm.
%
%
%
%
% Please run 'Example_iAPG_Step1.m' first.
%
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
close all

%% Factors after tensor decomposition
SpatialFactor	= A.U{1};
SpecTempFactor	= A.U{2};
SubjCondFactor	= A.U{3};

%% Label information
chanlocs = CurrentData.chanlocs; % Channel locations
ERP_frex = CurrentData.Fa; % Frequency
ERP_time = CurrentData.tim; % Time

%% Common files and paths
SplineFile='3D-64CHANNEL.SPL';

%% Image folder
CompImgPath=fullfile(pwd,'ResultImage');
if ~exist(CompImgPath,'dir'), mkdir(CompImgPath); end

%% Plot ERP components
for CompIndex = 1:R
    figure;
    set(gcf,'outerposition',get(0,'screensize'));% Full screen
    
    % Topography
    subplot(2,4,[1,2]);
    topoplot(zscore(SpatialFactor(:,CompIndex)),chanlocs);
    colorbar
    title(['#' int2str(CompIndex) ' Topography']);
    
    % 3D head topography
    subplot(2,4,[3,4]);
    headplot(zscore(SpatialFactor(:,CompIndex)),SplineFile,'electrodes','off');
    rotate3d on;
    title('3D Head Topography')
    
    % Spectrogram (Frequency-Time)
    Spec=reshape(SpecTempFactor(:,CompIndex),FrequencyMode,TimeMode);
    subplot(2,4,5);
    contourf(ERP_time,ERP_frex,Spec,100,'linecolor','none');
    set(gca,'YScale','log');
    set(gca,'YTick',logspace(log10(ERP_frex(1)),log10(ERP_frex(end)),15));
    set(gca,'YTickLabel',round(logspace(log10(ERP_frex(1)),log10(ERP_frex(end)),15)*10)/10);
    colorbar;
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    title('Spectrogram');
    
    % Subject-Condition
    subplot(2,4,7);
    bar(1:length(SubjCondFactor(:,CompIndex)),SubjCondFactor(:,CompIndex));
    xlim([0 length(SubjCondFactor(:,CompIndex))+1]);
    xlabel('Subject Index');
    ylabel('Strength');
    title('Subject-Condition');

    % Save figures
    sCompIndex=sprintf('%02d',CompIndex);
    saveas(gca,[CompImgPath filesep sCompIndex '.png'],'png');

end
