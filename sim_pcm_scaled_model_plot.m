%% Vizualize simulated "scaled" dataset from PCM
% Remi Gau - 2018-04-25
% Generates data using the PCM machinery (https://github.com/jdiedrichsen/pcm_toolbox)
% where condition 1 is a scaled version of condition 2 and then applies 
% a Z scoring row normalization to see the effect on the data.
% Results are vizualized to give an idea of how normalization might affect
% discriminability (by SVC on an MVPA analysis).
%
% The simulation are run with different values for the the 2 thetas of the
% scaled mode used to generate the data. Theta 2 is kept at 1, and theta 1
% increases from 1 to a set value in set amount of steps.

clc; clear;

StartDir = fullfile(pwd);
addpath(genpath(fullfile(StartDir, 'subfun')))
% adapt to point to wherever the PCM is on your machine
addpath('D:\Dropbox\GitHub\pcm_toolbox') 

Fig_dir = fullfile(StartDir, 'figures');
mkdir(Fig_dir)

numSim = 1; % number of subjects
NbVox = 500; % number of voxels (or vertices, channels...) in ROI
NbSess = 200; % number of fMRI sessions
MaxTheta1 = 6; % sets the upper limit of the theta range to try
NbSteps = 4; % number of steps between to try on the theta range
Model2Use = 1; % model to use to generate data (1 is the scaled model)


%% Define models
% In this case we only use the scaled model

% Scaled
Model{1}.type       = 'feature';
Model{end}.Ac = [1 0]';
Model{end}.Ac(:,1,2) = [0 1]';
Model{end}.name       = 'Scaled';
Model{end}.numGparams = size(Model{end}.Ac,3);
Model{end}.fitAlgorithm = 'NR';

% % Scaled and independent
% Model{end+1}.type       = 'feature';
% Model{end}.Ac = [1 0]';
% Model{end}.Ac(:,1,2) = [0 1]';
% Model{end}.Ac(:,2,3) = [0 1]';
% Model{end}.name       = 'Scaled+Independent';
% Model{end}.numGparams = size(Model{end}.Ac,3);
% Model{end}.fitAlgorithm = 'NR';
% 
% % Independent
% Model{end+1}.type       = 'feature';
% Model{end}.Ac = [1 0]';
% Model{end}.Ac(:,2,2) = [0 1]';
% Model{end}.name       = 'Independent';
% Model{end}.numGparams = size(Model{end}.Ac,3);
% Model{end}.fitAlgorithm = 'NR';


%% Set values to generate data with PCM machinery
%   theta:   numParams x 1 vector of parameters for Model
theta = [linspace(1,MaxTheta1,NbSteps)' ones(NbSteps,1)];

%   signal:  Signal variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox>
%   noise:   Noise  variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox>
% signal = ones(numSim,1)+randn(numSim,1);%rand(1,NbVox);
% noise = 1.5;%rand(1,NbVox);

signal = 1;%rand(1,NbVox);
noise = 1;%rand(1,NbVox);


%   numSim:  number of simulations,all returned in cell array Y

%   D: Experimental structure with fields
%       D.numPart = number of partititions
%       D.numVox  = number of independent voxels
D.numPart = NbSess;
D.numVox  = NbVox;

% VARARGIN:
%   'signalDist',fcnhnd:    Functionhandle to distribution function for signal
%   'noiseDist',fcnhnd:     Functionhandle to distribution function for noise
%   'design',X:             - Design matrix (for encoding-style models)
%                           - Condition vector (for RSA-style models)
%                           Design matrix and Condition vector are assumed
%                           to be for 1 partition only.
%                           If not specified - the function assumes a
%                           RSA-style model with G being numCond x numCond
noiseDist = @(x) norminv(x,0,1);   % Standard normal inverse for Noise generation
signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation

% Design matrix
X = [1 0;0 1];


%% Generate data and plot
close all

Nb_col_plot = 4;
FigDim = [50, 50, 1400, 700];
Visibility = 'on';

figure('name', 'PCM-MVPA', 'Position', FigDim, 'Color', [1 1 1], ...
    'Visible', Visibility);

iSubplot = 1;
for iTheta = 1:NbSteps
    
    %% Generate data
    [Y,partVec,condVec] = pcm_generateData(Model{Model2Use},theta(iTheta,:)',...
        D,numSim,signal,noise,'signalDist', noiseDist, ...
        'noiseDist', signalDist, 'design', X);

    G(:,:,iTheta) = pcm_estGCrossval(Y{1},partVec,condVec(:,1) + condVec(:,2)*2);
    
    %% Plots one condition against another for all voxels
    % (to show how one is the scaled version of the other)
    subplot(NbSteps,Nb_col_plot,iSubplot)
    hold on
    grid on
    
    % plot values of each voxel averaged over sessions
    scatter(...
        mean(Y{1}(logical(condVec(:,1)),:)),...
        mean(Y{1}(logical(condVec(:,2)),:)), 'b.')
    
    title(sprintf('Cdt 2 = f(Cdt 1)\nAveraged across sessions'))
    xlabel('Cdt 1')
    ylabel(sprintf('Theta 1 = %0.2f\n\nCdt 2',theta(iTheta,1)))
    ax = axis;
    axis([-20 20 ax(3) ax(4)])

    iSubplot = iSubplot + 1;
    
    
    %% plot one voxel against another for the 2 conditions across all sessions 
    %(to show "discriminability")
    subplot(NbSteps,Nb_col_plot,iSubplot)
    hold on
    grid on

    % Plot all examplars of condition 1 
    scatter(...
        Y{1}(logical(condVec(:,2)),1),...
        Y{1}(logical(condVec(:,2)),2), '.r')
    
    % Plot all examplars of condition 2 
    scatter(...
        Y{1}(logical(condVec(:,1)),1),...
        Y{1}(logical(condVec(:,1)),2), '.b')    
    
    title(sprintf('Voxel 1 = f(Voxel 2)\nAll sessions'))
    xlabel('Voxel 1')
    ylabel('Voxel 2')
    
    if iTheta == 1
    legend({'Cdt 1','Cdt 2'},'Location','SouthEastOutside')
    end

    iSubplot = iSubplot + 1;
    

    %% same as above but after Z scoring across each row (examplar)
    tmp = zscore(Y{1},0,2);
    G_Zscore(:,:,iTheta) = pcm_estGCrossval(tmp,partVec,condVec(:,1) + condVec(:,2)*2);
    
    %%
    subplot(NbSteps,Nb_col_plot,iSubplot)
    hold on
    grid on

    scatter(...
        mean(tmp(logical(condVec(:,1)),:)),...
        mean(tmp(logical(condVec(:,2)),:)), 'b.')
    
    title(sprintf('Cdt 2 = f(Cdt 1)\nrow Z scoring'))
    xlabel('Cdt 1')
    ylabel('Cdt 2')

    iSubplot = iSubplot + 1;
    
    %% 
    subplot(NbSteps,Nb_col_plot,iSubplot)
    hold on
    grid on

    scatter(...
        tmp(logical(condVec(:,2)),1),...
        tmp(logical(condVec(:,2)),2), '.r')
    scatter(...
        tmp(logical(condVec(:,1)),1),...
        tmp(logical(condVec(:,1)),2), '.b')    
    
    title(sprintf('Voxel 1 = f(Voxel 2)\nrow Z scoring'))
    xlabel('Voxel 1')
    ylabel('Voxel 2')

    iSubplot = iSubplot + 1;
    

end

print(gcf, fullfile(Fig_dir, 'PCM_Scaled.tif'), '-dtiff')


%% Plot the corresponding G matrices
close all

FigDim = [50, 50, 600, 700];

figure('name', 'PCM-MVPA-GMat', 'Position', FigDim, 'Color', [1 1 1], ...
    'Visible', Visibility);

CLIM = [0 max(G(:))];
CLIM2 = [0 max(G_Zscore(:))];

iSubplot = 1;
for iTheta = 1:NbSteps
    subplot(NbSteps,2,iSubplot)
    imagesc(G(:,:,iTheta),CLIM)
    set(gca,'xtick', 1:2, 'xticklabel', ({'Cdt 1', 'Cdt 2'}),...
        'ytick', 1:2, 'yticklabel', ({'Cdt 1', 'Cdt 2'}) );
    axis square
    title(sprintf('theta 1 = %3.2f ; theta 2 = 1', theta(iTheta,1)))
    iSubplot = iSubplot + 1;
    
    subplot(NbSteps,2,iSubplot)
    imagesc(G_Zscore(:,:,iTheta),CLIM2)
    set(gca,'xtick', 1:2, 'xticklabel', ({'Cdt 1', 'Cdt 2'}),...
        'ytick', 1:2, 'yticklabel', ({'Cdt 1', 'Cdt 2'}) );
    axis square
    title(sprintf('theta 1 = %3.2f ; theta 2 = 1\n Z-score normalization', theta(iTheta,1)))
    iSubplot = iSubplot + 1;
end

print(gcf, fullfile(Fig_dir, 'PCM_Scaled_G_Mat.tif'), '-dtiff')