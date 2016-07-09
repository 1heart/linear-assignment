%--------------------------------------------------------------------------
% Function:    distMatsAcrossLevels
% Description: Creates all the necessary Euclidean distance matrices for 
%              the different wavelet decomposition levels.  
%
% Inputs:
%
% Outputs:
%   distMats          - Lx1 cell array where L represents the
%                       number of decomposition levels. distMat{i} returns
%                       an (N_i x N_i) matrix, where N_i is the number 
%                       of distance locations evaluated at the ith level.
%
% Usage:
%
% Authors(s):
%   Adrian M. Peter
%--------------------------------------------------------------------------
function [distMats basesLocs]  = euclidDistsBetwCoeffs3D(wName, startLevel,...
                                 stopLevel, coeffsIdx, sampleDataSupp)

distMats  = cell(length(startLevel:stopLevel),1);
basesLocs = [];

% Loop over all levels and compute distance matrices for each level.
l = 0;
for j = startLevel : stopLevel
    l        = l + 1;
    transX   = translationRange(sampleDataSupp(1,:), wName, j);
    transY   = translationRange(sampleDataSupp(2,:), wName, j);
    transZ   = translationRange(sampleDataSupp(3,:), wName, j);
    wSupport = waveSupport(wName);
    
    % Get end point of basis support for the current level.
    basisEndPt = wSupport(2)/2^startLevel;
    
    % Get the centroid position of the basis for the current level.
    basisCenter = basisEndPt/2;
    
    % Now get the position of all the basis for all the translation values.
    basisLocsX = [transX(1):transX(2)]'*(1/2^startLevel)+basisCenter;
    basisLocsY = [transY(1):transY(2)]'*(1/2^startLevel)+basisCenter;
    basisLocsZ = [transZ(1):transZ(2)]'*(1/2^startLevel)+basisCenter;
    
    if(l==1)
        numCoeffs = diff(coeffsIdx(l,:))+1;
    else
        numCoeffs = diff(coeffsIdx(l+1,:))+1;
        numCoeffs = numCoeffs/3;
    end
    if( (length(basisLocsX)*length(basisLocsY))~=numCoeffs )
        error('Coefficient length does not match number of basis location');
    end
    
    % Matrix of all basis locations.
    [xx, yy, zz] = meshgrid(basisLocsX, basisLocsY, basisLocsZ);
    basisLocs   = [xx(:), yy(:), zz(:)];
    numLocs     = size(basisLocs,1);
    % Create distance matrix.
    cVec1       = repmat(basisLocs,numLocs,1);
    temp        = repmat(basisLocs',numLocs,1);
    cVec2       = reshape(temp,2,numLocs^2)';
    % Trying out square of the distance.
    distMats{l} = reshape((sum((cVec1-cVec2).^2,2)),numLocs,numLocs);
    
    tempDistMat = pdist2(cVec1,cVec2);
    
    norm(tempDistMat - distMats{1}, 'fro')
    
    %distMats{l} = reshape(sqrt(sum((cVec1-cVec2).^2,2)),numLocs,numLocs);
    basesLocs   = [basesLocs basisLocs'];
end
basesLocs = basesLocs';
