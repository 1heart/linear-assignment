%--------------------------------------------------------------------------
% Function:    linearAssignment
% Description: Performs linear assignment between wavelet coefficients from
%              shape 1 versus shape 2.  Calls C++ routines to do actual
%              work.  The C++ codes implements:
%                A Shortest Path Algorithm for Dense and Sparse Linear 
%                Assignment Problems (R. Jonker), Computing 38 
%                (1987) 325-340.
%              The code is available at:
%                http://www.magiclogic.com/assignment.html
%
% Inputs:
%   costMatrix        - NxN matrix of costs for each translation location
%                       over the support of the shape.
%   weights           - NxN matrix of penalty weights for each grid
%                       location where the cost is defined.
%   namePostfix       - any character string that needs to put after each
%                       shape file name.
%   lambda            - importance weighting for Eulcidean distance matrix
%                       used to penalize large swaps during linear
%                       assignment. If lambda=0 there will not be a
%                       penalty. 
%   addWeights        - Flag to indicate if we want to include the weights
%                       in the cost.  Could accomplish the same effect with
%                       lambda=0, but has already been written ;).
%   maximize          - Used to indicate if we want to maximize or minimize
%                       the cost.  Necessary there that can be a global
%                       sign difference between shapes.
%   rewriteFile       - Used to indicate if we need to rewrite cost matrix
%                       text file.
% Outputs:
%   assignments       - vector containing the new coefficient ordering.
%
% Usage:
%
% Authors(s):
%   Adrian M. Peter and Anand Rangarajan
%--------------------------------------------------------------------------
function assignments = linearAssignment(costMatrix, weights, lambda,...
                                        addWeights, maximize, ...
                                        inputFileName, outputFileName,...
                                        rewriteFile, execPath)
if(~exist('inputFileName','var'))
    inputFileName  = 'costMatrix.cm';
end

if(~exist('outputFileName','var'))
    outputFileName = 'assignments.txt';
end

if(~exist('rewriteFile','var'))
    % Flag to indicate if we need to rewrite cost matrix text file.
    rewriteFile = 1;
end

if(maximize)
    costMatrix = -costMatrix;
end

% Scale, shift, round matrix cost matrix in order to get linear assignment
% code to converge.  We are adding the weights b/c the LAP code does a
% minimization instead of maximizing.  
costMatrix = costMatrix+lambda*(weights);
% Right now only integer value costs seem to conver all the time.  Just
% have to make sure that our cost matrix values do not exceed the maximum
% integer values.  On MS platform:
%       Minimum value for a variable of type int. 	�2147483647 � 1
%       Maximum value for a variable of type int. 	2147483647 [2^(31)-1]
costMatrix = round(costMatrix);
if(min(costMatrix(:))<-2147483648)
    error('Minimum cost values exceed machine integer precision');
end
if(max(costMatrix(:))>2147483647)
    error('Maximum cost values exceed machine integer precision');
end

%posMat = posMat - min(posMat(:));

% FIX: May have to subtract weights if we are minimizing.
if(rewriteFile)
%    if(addWeights)
 %       writeMatrixToFile(posMat,inputFileName);
        %writeMatrixToFile(costMatrix+lambda*(weights),inputFileName);
%    else
%        writeMatrixToFile(costMatrix,inputFileName);
        writeMatrixToBinaryFile(costMatrix,inputFileName);
%    end
end

%[status shellText] = dos(['..\lap_cpp\Release\linAssign.exe .\' inputFileName ' .\' outputFileName]);
if(ispc)
%     [status shellText] = dos(['.\lap_binInFile\Release\linAssignBinFile.exe .\' inputFileName ' .\' outputFileName]);
%     [status shellText] = dos(['C:\Users\mmoyou\Documents\MATLAB\CurrentResearch\Wasserstein\ShapeLaneRouge\SlidingSqWDEDist\lap_binInFile\Release\linAssignBinFile.exe '...
%                         inputFileName ' ' outputFileName]);
    [status shellText] = dos([execPath, 'linAssignBinFile.exe '...
                        inputFileName ' ' outputFileName]);

else
    [status shellText] = unix(['./linAssignBinFile ' inputFileName ' ./' outputFileName]);
end

% try
    load(['./' outputFileName]);
% catch
%     a = 1;
% end
% noExtName = outputFileName(1:findstr(outputFileName,'.')-1);
% noExtName = outputFileName(1:findstr(outputFileName,'.txt')-1);



initNameInd = strfind(outputFileName, 'assignments');
endNameInd = strfind(outputFileName, '.txt');
evalAssignName = outputFileName(initNameInd : (endNameInd - 1));

eval(['assignments = ' evalAssignName ';']);
delete(inputFileName);
delete(outputFileName);

% eval(['assignments = ' noExtName ';']);
% eval(['assignments = ' noExtName(15:end) ';']);

%load assignments.txt;
% assignments
% e2(assignments+1) 
