%% Simulated "scaled" dataset for PCM
% Remi Gau - 2018-04-25

clc; clear; close all

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox')

numSim = 1; % number of subjects
NbVox = 500; % number of voxels or vertices in ROI
NbSess = 200; % number of fMRI sessions
theta1 = 6; % sets the upper limit of the theta range to try
NbSteps = 4; % number of steps between to try on the theta range
Model2Use = 1; % model to use to generate data


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
theta = [linspace(1,theta1,NbSteps)' ones(NbSteps,1)];

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
FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

figure('name', 'PCM-MVPA', 'Position', FigDim, 'Color', [1 1 1], ...
    'Visible', Visibility);

iSubplot = 1;
for iTheta = 1:NbSteps
    
    %% Generate data
    [Y,partVec,condVec] = pcm_generateData(Model{Model2Use},theta(iTheta,:)',...
        D,numSim,signal,noise,'signalDist', noiseDist, ...
        'noiseDist', signalDist, 'design', X);
    
    %% Plots one condition against another for all voxels
    % (to show how one is the scaled version of the other)
    subplot(NbSteps,Nb_col_plot,iSubplot)
    hold on
    grid on
    
    % plot values of each voxel averaged over sessions
    scatter(...
        mean(Y{1}(logical(condVec(:,1)),:)),...
        mean(Y{1}(logical(condVec(:,2)),:)), 'b.')
    
    title('Cdt 2 = f(Cdt 1)')
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
    
    title('Voxel 1 = f(Voxel 2)')
    xlabel('Voxel 1')
    ylabel('Voxel 2')

    iSubplot = iSubplot + 1;
    
    
    
    
    %% same as above but after Z scoring across each row (examplar)
    tmp = zscore(Y{1},0,2);
    
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

