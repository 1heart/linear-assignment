% Testing shape Lane Rouge.
clc; close all; clear;

addpath('C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\WaveletDensityEstimator\wde2D\');
addpath('C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\WaveletDensityEstimator\waveCommon\');
addpath('.\lap_binInFile\');
addpath('C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\Hypersphere_Code');

% Load the coefficients. 
load('testCoeffs2');

saveFold = 'C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\ShapeLaneRouge\SlidingSqWDEDist\ShapeLaneRougeDistMatrices\Brown\';

% % Compute the Karcher mean of th 
% meanMat = computeKarcherMeanForDataset(11,9, coeffs);

% Get the number of coefficients. 
numCoeffs = size(coeffs,1);

slrParamLb = 1;
slrParamUb = 2000;
numSlrVals = 10;
slrVals = linspace(slrParamLb, slrParamUb, numSlrVals);
% slrParam = 2000;
slrDistCell = cell(numSlrVals,1);

F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F);
h = waitbar(0,'Please wait...');
totalComps = numCoeffs*numCoeffs;
k = 0;

coeffsIdx = coefficientIndex(wdeSet.startLevel, wdeSet.stopLevel,...
                                      wdeSet.sampleSupp, wdeSet.wName,...
                                      1);
distMatInd = reshape(1:numCoeffs^2, [numCoeffs, numCoeffs]);
execLoc = 'C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\ShapeLaneRouge\SlidingSqWDEDist\lap_binInFile\Release\';

% Create a matrix of distances between the coefficients. 
for slrInd = 1 : numSlrVals
    slrDist = zeros(numCoeffs);
    slrParam = slrVals(slrInd);
    for i = 1 : numCoeffs
        i
        coeffsi = coeffs{i,1};
        st = tic;
        for j = 1 : numCoeffs
            if (j > i)
%                 j
                slrDist(i,j) = wde2DSlidingDist(coeffsi, coeffs{j,1}, wdeSet.wName,...
                 wdeSet.startLevel,wdeSet.stopLevel, coeffsIdx,...
                     wdeSet.sampleSupp,slrParam,...
                     ['.\CostMatrices\costMatrix', num2str(j), '.cm'],...
                     ['.\Assignments\assignments', num2str(j), '.txt'],...
                     execLoc);
            
            
%                 slrDist(i,j) = wde2DSlidingDist(coeffsi, coeffs{j,1}, wdeSet.wName,...
%                  wdeSet.startLevel,wdeSet.stopLevel, coeffsIdx,...
%                      wdeSet.sampleSupp,slrParam,...
%                      ['C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\ShapeLaneRouge\SlidingSqWDEDist\CostMatrices\CostMatrices\costMatrix', num2str(j), '.cm'],...
%                      ['C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\ShapeLaneRouge\SlidingSqWDEDist\Assignments\Assignments\assignments', num2str(j), '.txt']);
            
            end
        end
        stopTime = toc(st);
        disp(['Stop time = ', num2str(stopTime)]);
    end
    
    slrDistCell{slrInd,1} = slrDist + slrDist';
end
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F);

saveName = [saveFold, 'brown_' num2str(slrParam)];
save(saveName, 'slrDistCell');

%     
% 
% c1Num = 1;
% c2Num = 20;
% 
% c1 = coeffs{c1Num,1};
% c2 = coeffs{c2Num,1};
% 
% % Compute the coeffIdx.
% coeffsIdx = coefficientIndex(wdeSet.startLevel, wdeSet.stopLevel,...
%                                       wdeSet.sampleSupp, wdeSet.wName,...
%                                       1);
% 
% % Call the sliding wavelet distance function.
% slidingDist = wde2DSlidingDist(c1, c2, wdeSet.wName,...
%                      wdeSet.startLevel,wdeSet.stopLevel, coeffsIdx,...
%                      wdeSet.sampleSupp,2000)
% 
% distMat = distBtwCoefficients(c1,c2)                 