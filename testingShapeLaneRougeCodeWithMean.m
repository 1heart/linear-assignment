% Testing shape Lane Rouge.
clc; close all; clear;

% warning off;
% addpath('C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\WaveletDensityEstimator\wde2D\');
% addpath('C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\WaveletDensityEstimator\waveCommon\');
% addpath('.\lap_binInFile\');
% addpath('C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\Hypersphere_Code');

% Load the coefficients. 
load('testCoeffs2_v3');

execLoc = 'C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\ShapeLaneRouge\SlidingSqWDEDist\lap_binInFile\Release\';

% Get the number of coefficients. 
[numCoeffs, numShapes] = size(meanMat);

slrParamLb = 10;
slrParamUb = 3000;
numSlrVals = 10;
slrVals = linspace(slrParamLb, slrParamUb, numSlrVals);

numShapesm1 = numShapes - 1; % Number of shapes minus 1.
% Store the distance for each slr parameter. 
slrDistMat = zeros(numSlrVals, numShapesm1); 


coeffsIdx = coefficientIndex(wdeSet.startLevel, wdeSet.stopLevel,...
                                      wdeSet.sampleSupp, wdeSet.wName,...
                                      1);
shapeNum = 2;                                  
% Distance between the first coefficient and the rest.
coeffDists = distBtwCoefficients(coeffMat(:,shapeNum), meanMat);                                 
% Create a matrix of distances between the coefficients. 
for slrInd = 1 : numSlrVals
    slrParam = slrVals(slrInd); % Current slr parameter. 
    st = tic;
    parfor j = 2 : numShapes
         slrDistMat(slrInd,j-1) = wde2DSlidingDist(coeffMat(:,shapeNum),...
                        meanMat(:,j),...
                        wdeSet.wName,wdeSet.startLevel,wdeSet.stopLevel,... 
                        coeffsIdx, wdeSet.sampleSupp, slrParam,...
             ['.\CostMatrices\costMatrix', num2str(j), '.cm'],...
             ['.\Assignments\assignments', num2str(j), '.txt'],execLoc);

    end
    disp(slrParam);
    disp(num2str(coeffDists));
    disp(num2str(slrDistMat(slrInd,:)));
    stopTime = toc(st);
    disp(['Stop time = ', num2str(stopTime)]);
end
               