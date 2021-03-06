%% Run MVPA on simulated "scaled" dataset from PCM
% Remi Gau - 2018-04-25
% Generates data using the PCM machinery (https://github.com/jdiedrichsen/pcm_toolbox)
% where condition 1 is a scaled
% version of condition 2 and then runs an SVC on it with different type of
% normalization (e.g Z scoring row normalization and/or mean centering
% column normalization)
%
% The simulation are run with different values for the the 2 thetas of the
% scaled mode used to generate the data. Theta 2 is kept at 1, and theta 1
% increases from 1 to a set value in set amount of steps.
% 
% Uses the libsvm for the SVC. Might require recompiling of some of the mex
% function if you are not running this on windows (c and make files are in
% subfun/libsvm/matlab)

clc; clear; close all

StartDir = fullfile(pwd);
addpath(genpath(fullfile(StartDir, 'subfun')))
% adapt to point to wherever the PCM is on your machine
addpath('D:\Dropbox\GitHub\pcm_toolbox') 

Save_dir = fullfile(StartDir, 'results');
mkdir(Save_dir)
Fig_dir = fullfile(StartDir, 'figures');
mkdir(Fig_dir)


NbSim = 10; % number of simulations (or subjects)
NbVox = 200; % number of voxels / channels / features (in MVPA lingo)
NbSess = 20; % number of examplars for each condition

MaxTheta1 = 6; % max value taken by theta 1
NbSteps = 5; % number of steps between theta1 = 1 and its max value



%% Define models
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


%% Set values for PCM data generation
%   theta:   numParams x 1 vector of parameters for Model
theta = [linspace(1,MaxTheta1,NbSteps)' ones(NbSteps,1)];

%   signal:  Signal variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox>
%   noise:   Noise  variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox>
signal = ones(NbSim,1)+randn(NbSim,1);%rand(1,NbVox);
noise = 1.5;%rand(1,NbVox);

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


%% Set MVPA parameters

% Feature selection (FS)
opt.fs.threshold = 0.75;
opt.fs.type = 'ttest2';

% Recursive feature elminiation (RFE)
opt.rfe.threshold = 0.01;
opt.rfe.nreps = 20;

% SVM C/nu parameters and default arguments
opt.svm.machine = 'C-SVC';
opt.svm.log2c = 1;
opt.svm.dargs = '-s 0';

opt.svm.dargs = [opt.svm.dargs ' -t 0 -q']; % inherent linear kernel, quiet mode

opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.scaling.idpdt = 1; % scale test and training sets independently
opt.permutation.test = 0; % label permutation 

% Choices of image (row) or voxel/feature (column) normalization
opt.scaling.img.eucledian = 0;
opt.scaling.feat.mean = 0;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;

% Necessary to tell the MVPA functions what SVC to run
SVM(1) = struct('name', 'Cst 1 VS Cdt 2', 'class', [1 2], 'ROI_2_analyse', 1);

% CV scheme (based on leaving one run out from 3 different days)
CV_id = 1:NbSess;
sets = {1:6,1:7,1:7};
[x, y, z] = ndgrid(sets{:});
TestSessList{1,1} = [x(:) y(:) z(:)];
NbCV = size(TestSessList{1,1}, 1);


%% START

for iTheta = 1:size(theta,1)
    
    % generate data
    [Y,partVec,condVec] = pcm_generateData(Model{1},theta(iTheta,:)',D,NbSim,signal,noise,...
        'signalDist', noiseDist, 'noiseDist', signalDist, 'design', X);
    
    % to organize cross validation for MVPA
    CV_Mat(:,1) = condVec(:,1) + condVec(:,2)*2;
    CV_Mat(:,2) = partVec;
    
    %%
    for iSim = 1:NbSim
        
        for iCV=1:NbCV
            
            TestSess = []; %#ok<NASGU>
            TrainSess = []; %#ok<NASGU>
            
            % Separate training and test sessions
            [TestSess, TrainSess] = deal(false(size(1:NbSess)));
            
            TestSess(TestSessList{1,1}(iCV,:)) = 1; %#ok<*PFBNS>
            TrainSess(setdiff(CV_id, TestSessList{1,1}(iCV,:)) )= 1;
            
            % Run MVPA with no normalization
            opt.scaling.img.zscore = 0;
            opt.scaling.feat.mean = 0;
            results = machine_SVC(SVM(1), Y{iSim}, CV_Mat, TrainSess, TestSess, opt);
            Acc(iSim,1,iCV,iTheta) = mean(results.pred==results.label);
            
            % Run MVPA with Z scoring for row normalization (across voxels)
            opt.scaling.img.zscore = 1;
            opt.scaling.feat.mean = 0;
            results = machine_SVC(SVM(1), Y{iSim}, CV_Mat, TrainSess, TestSess, opt);
            Acc(iSim,2,iCV,iTheta) = mean(results.pred==results.label);
            
            % Run MVPA with mean centering for row normalization (across
            % examplars)
            opt.scaling.img.zscore = 0;
            opt.scaling.feat.mean = 1;
            results = machine_SVC(SVM(1), Y{iSim}, CV_Mat, TrainSess, TestSess, opt);
            Acc(iSim,3,iCV,iTheta) = mean(results.pred==results.label);
            
            % Run MVPA with both normalization ()
            opt.scaling.img.zscore = 1;
            opt.scaling.feat.mean = 1;
            results = machine_SVC(SVM(1), Y{iSim}, CV_Mat, TrainSess, TestSess, opt);
            Acc(iSim,4,iCV,iTheta) = mean(results.pred==results.label);
            
        end
    end
    
end

%% Saves data
clear iCV iSubj iTheta Y

save(fullfile(Save_dir, ['PCM_MVPA_', datestr(now, 'yyyy_mm_dd_HH_MM'), '.mat']))
