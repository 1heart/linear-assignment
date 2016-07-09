%--------------------------------------------------------------------------
% Function:    wde2DSlidingDist
% Description: Computes the cosine distance between two square-root wavelet 
%              densities (sqWDE) after "adjusting" one coefficient set to 
%              be more similar to the other.
%
% Inputs:
%   coeffs1           - Coefficient vector of first sqWDE.
%   coeffs2           - Coefficient vector of first sqWDE.
%   wName             - Name of wavelet used during density estimation.
%                       Use matlab naming convention for wavelets.
%   startLevel        - Starting level for the the father wavelet
%                       (i.e. scaling function).  
%   stopLevel         - Last level for mother wavelet scaling.  The start
%                       level is same as the father wavelet's.
%   coeffsIdx         - Lx2 matrix containing the start and stop index
%                       locations of the coeffients for each level in the
%                       coefficient vector.  L is the number of levels.
%   sampleDataSupport - 2x2 matrix of the sample support.
%                       First row gives min x value and max x value
%                       Second row gives min y value and max y value
%   lambda            - importance weighting for Eulcidean distance matrix
%                       used to penalize large swaps during linear
%                       assignment. If lambda=0 there will not be a
%                       penalty. (optional)
% Outputs:
%   slidingDist       - Dist on unit hypersphere after moving coefficients.
%
% Usage:
%
% Authors(s):
%   Adrian M. Peter
%--------------------------------------------------------------------------
function slidingDist = wde2DSlidingDist_Distributed(coeffs1, coeffs2,...
                            wName, startLevel, ...
                            stopLevel, coeffsIdx, sampleDataSupp,lambda,...
                            inputFileNameCostMat, outputFileNameAssign,...
                            execLoc)

if(~exist('lambda'))
    lambda = 5e-4;
end

if(lambda == 0)    
    addWeights = 0; % Flag to signal use of weights.
else
    addWeights = 1;
end

perturbCosts = 0; % Flag to perturb cost matrix to break symmetry. 
                  % (FIX: May have to just scale.)      

% Call function to get cell array of all the necessary Euclidean distance
% matrices for the different levels.  We will have to load only one of the
% shape mat files since all of them are required to have the same level of
% decompositions.
[distMats, ~] = euclidDistsBetwCoeffs3D(wName, startLevel,...
                                 stopLevel, coeffsIdx, sampleDataSupp);

numLevels = (stopLevel-startLevel)+1;

if( (numLevels == 1) && (size(coeffsIdx,1)==1) )
    scalingOnly = 1;
else
    scalingOnly = 0;
end

% Go through all levels for the current shape pair and solve
% assignment problem to get correct coefficient ordering.
shape1           = coeffs1;
shape2           = coeffs2;
saveCoeffsOrder1 = [];
saveCoeffsOrder2 = [];
for j = 1 : numLevels
    costMatrix = [];
    % Need to get correct coefficient locations for wavelet basis
    % since there are 3 per x-y location.
    if(~scalingOnly)
        numCoeffsPerBasis = (coeffsIdx(j+1,2)-coeffsIdx(j+1,1)+1)/3;
        b1CoeffsIdxStart  = coeffsIdx(j+1,1);
        b1CoeffsIdxEnd    = coeffsIdx(j+1,1)+(numCoeffsPerBasis-1);
        b2CoeffsIdxStart  = b1CoeffsIdxEnd + 1;
        b2CoeffsIdxEnd    = b2CoeffsIdxStart+(numCoeffsPerBasis-1);
        b3CoeffsIdxStart  = b2CoeffsIdxEnd + 1;
        b3CoeffsIdxEnd    = b3CoeffsIdxStart+(numCoeffsPerBasis-1);
        if(b3CoeffsIdxEnd ~= coeffsIdx(j+1,2))
            error(['Error calculating wavelet coefficient indices at level: '...
                   num2str(j+1)]);
        end
    end

    % This first level can have scaling and wavelet bases.
    if(j==1)
        s1ScalingCoeffs = shape1(coeffsIdx(j,1):coeffsIdx(j,2));
        s2ScalingCoeffs = shape2(coeffsIdx(j,1):coeffsIdx(j,2));
        % Create cost matrix as outer product of wavelet coefficients.
        costMatrix = s1ScalingCoeffs*s2ScalingCoeffs';
        % If wavelet coeffs are present need to add them to cost.
        if(~scalingOnly)
            % Need to add 3 cost matrices together b/c we have 3
            % wavelet coefficient per x-y grid location.
            s1WaveletCoeffs1 = shape1(b1CoeffsIdxStart:b1CoeffsIdxEnd);
            s2WaveletCoeffs1 = shape2(b1CoeffsIdxStart:b1CoeffsIdxEnd);
            costMatrix       = costMatrix + s1WaveletCoeffs1*s2WaveletCoeffs1';

            s1WaveletCoeffs2 = shape1(b2CoeffsIdxStart:b2CoeffsIdxEnd);
            s2WaveletCoeffs2 = shape2(b2CoeffsIdxStart:b2CoeffsIdxEnd);
            costMatrix       = costMatrix + s1WaveletCoeffs2*s2WaveletCoeffs2';

            s1WaveletCoeffs3 = shape1(b3CoeffsIdxStart:b3CoeffsIdxEnd);
            s2WaveletCoeffs3 = shape2(b3CoeffsIdxStart:b3CoeffsIdxEnd);
            costMatrix       = costMatrix + s1WaveletCoeffs3*s2WaveletCoeffs3';
        end
    else % j>1
        s1WaveletCoeffs1 = shape1(b1CoeffsIdxStart:b1CoeffsIdxEnd);
        s2WaveletCoeffs1 = shape2(b1CoeffsIdxStart:b1CoeffsIdxEnd);
        costMatrix       = s1WaveletCoeffs1*s2WaveletCoeffs1';

        s1WaveletCoeffs2 = shape1(b2CoeffsIdxStart:b2CoeffsIdxEnd);
        s2WaveletCoeffs2 = shape2(b2CoeffsIdxStart:b2CoeffsIdxEnd);
        costMatrix       = costMatrix + s1WaveletCoeffs2*s2WaveletCoeffs2';

        s1WaveletCoeffs3 = shape1(b3CoeffsIdxStart:b3CoeffsIdxEnd);
        s2WaveletCoeffs3 = shape2(b3CoeffsIdxStart:b3CoeffsIdxEnd);
        costMatrix       = costMatrix + s1WaveletCoeffs3*s2WaveletCoeffs3';
    end % if(j==1)
    
    % Compute linear assignment.  We have to do it twice and see
    % which one maximizes the inner prodcut between the shape
    % coefficients.  (Shapes can have global sign difference
    % between them since they were estimated as sqrt of the
    % density.
    costScale = 10000/max([(shape1.^2)' (shape2.^2)']);
    disp(['Scaling costmatrix by: ' num2str(costScale)]);
    costMatrix  = costScale*costMatrix;
    
    % We will compute the + and - version of the cost to determine which
    % gives us the best ordering.
    maximize    = 1;
    t0          = clock; 
%     coeffOrder1 = linearAssignment(costMatrix, distMats{j}, lambda, addWeights, maximize);
    coeffOrder1 = linearAssignment_Distributed(...
                costMatrix, distMats{j}, lambda, addWeights, maximize,...
                  inputFileNameCostMat, outputFileNameAssign,1,execLoc);
    sec         = etime(clock,t0);
    if(sec <= 60)
        disp(['LA1 took: ' num2str(sec) ' seconds.'])
    elseif((sec > 60) && (sec <= 3600))
        disp(['LA1 took: ' num2str(sec/60) ' minutes.'])
    else
        disp(['LA1 took: ' num2str(sec/3600) ' hours.'])
    end

    if(length(coeffOrder1)~=size(costMatrix,1))
        error(['LA1 did not return proper ordering vector!']);
    end

    maximize    = 0;
    t0          = clock; 
    coeffOrder2 = linearAssignment_Distributed(...
                costMatrix, distMats{j}, lambda, addWeights, maximize,...
                  inputFileNameCostMat, outputFileNameAssign,1,execLoc);
    sec         = etime(clock,t0);
    if(sec <= 60)
        disp(['LA2 took: ' num2str(sec) ' seconds.'])
    elseif((sec > 60) && (sec <= 3600))
        disp(['LA2 took: ' num2str(sec/60) ' minutes.'])
    else
        disp(['LA2 took: ' num2str(sec/3600) ' hours.'])
    end
    if(length(coeffOrder2)~=size(costMatrix,1))
        error(['LA2 did not return proper ordering vector!']);
    end

    % Save the coefficient ordering.
    if(j==1)
        % When saving the wavelet coefficient ordering we need to
        % add the correct offset so all 3 basis per x-y location are
        % reordered properly.
        if(~scalingOnly)
            saveCoeffsOrder1 = [saveCoeffsOrder1 (coeffOrder1+coeffsIdx(j,1))'...
                                (coeffOrder1+b1CoeffsIdxStart)'...
                                (coeffOrder1+b2CoeffsIdxStart)'...
                                (coeffOrder1+b3CoeffsIdxStart)']; 
            saveCoeffsOrder2 = [saveCoeffsOrder2 (coeffOrder2+coeffsIdx(j,1))'...
                                (coeffOrder2+b1CoeffsIdxStart)'...
                                (coeffOrder2+b2CoeffsIdxStart)'...
                                (coeffOrder2+b3CoeffsIdxStart)']; 
        else
            saveCoeffsOrder1 = [saveCoeffsOrder1 (coeffOrder1+coeffsIdx(j,1))'];
            saveCoeffsOrder2 = [saveCoeffsOrder2 (coeffOrder2+coeffsIdx(j,1))'];
        end
    else
        saveCoeffsOrder1 = [saveCoeffsOrder1 (coeffOrder1+b1CoeffsIdxStart)'...
                            (coeffOrder1+b2CoeffsIdxStart)'...
                            (coeffOrder1+b3CoeffsIdxStart)']; 
        saveCoeffsOrder2 = [saveCoeffsOrder2 (coeffOrder2+b1CoeffsIdxStart)'...
                            (coeffOrder2+b2CoeffsIdxStart)'...
                            (coeffOrder2+b3CoeffsIdxStart)']; 
    end

end % for j = 1 : numLevels

% Get two candidate shape orderings.  Remember we have two
% candidate orderings due to the possibility of have a global sign
% difference between the shapes.
shape2C1 = shape2(saveCoeffsOrder1);
shape2C2 = shape2(saveCoeffsOrder2);

% Get inner product for both candidates and pick the greater one.
innerProd1LA = shape1'*shape2C1;
innerProd2LA = shape1'*shape2C2;

% Compute shape distance on the sphere.
if(innerProd1LA >= innerProd2LA)
    slidingDist = acos(innerProd1LA);
else
    warining(['Inner product 2 with LA higher, ' num2str(innerProd2LA)]);
    slidingDist = acos(innerProd2LA);
end

if(~isreal(slidingDist))
    warning(['Sliding distances is complex with, Real = ' num2str(real(slidingDist)) ' Img = ' num2str(imag(slidingDist)) '. Returning real part.']);
    slidingDist = real(slidingDist);
end    
