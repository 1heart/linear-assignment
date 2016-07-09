function distMat = computeShapeDistance2D(X, Y, distType, varargin)
% distMat = computeShapeDistance2D(X, Y, distType)
% X is some DxN data.
% Y is some DxM data
X = real(X);
Y = real(Y);

switch distType
    case 'acos'
        distMat = real(acos(real(X)' * real(Y)));
    case 'slr'
        
        if (isempty(varargin) == 1)
            error('Provide the wavelet type data structure');
        else       
            wdeSet = varargin{1};
        end
        
        % Wavelet parameters. 
        wName = wdeSet.wName;
        startLevel = wdeSet.startLevel;
        stopLevel = wdeSet.stopLevel; 
        sampleDataSupp = wdeSet.sampleSupp;
        slrParam = wdeSet.slrParam;
        costMatFoldPath = wdeSet.costMatFoldPath;
        assignMatFoldPath = wdeSet.assignMatFoldPath;
        execPath = wdeSet.execPath;
        
        coeffsIdx = coefficientIndex(wdeSet.startLevel, wdeSet.stopLevel,...
                                      wdeSet.sampleSupp, wdeSet.wName,...
                                      1);
        
        % Check if X is equal to Y then the distance matrix will be
        % symmetric and you only have to compute the upper triangular half.
        isXEqualY = isequal(X,Y);
        
        numTestSamps = size(X,2); % Number of test samples. 
        numTrainSamps = size(Y,2); % Number of train samples. 
                
        % Initialize the distance matrix. 
        distMat = zeros(numTestSamps, numTrainSamps); 
        
        h = waitbar(0,'Computing...');
        % Loop through the shapes and compute the entries to the distance
        % matrix. 
        if (isXEqualY == 1) % Distance matrix is symmetric.
            for i = 1 : numTestSamps
                waitbar(i/numTestSamps,h)
                currX = X(:,i); % Get the current coefficient vector.
                parfor k = 1 : numTrainSamps
                    if (k > i)
                        k
                        % Compute the i,j distance. 
                        distMat(i,k) = wde2DSlidingDist(...
                                 currX, Y(:,k),...
                                 wName, startLevel, ...
                                 stopLevel, coeffsIdx, sampleDataSupp,...
                                 slrParam,...
                        [costMatFoldPath, 'costMatrix', num2str(k), '.cm'],...
                        [assignMatFoldPath, 'assignments', num2str(k), '.txt'],...
                        execPath);
                    end
                end
            end
            
            % Create a symmetric matrix with zeros on the diagonal.
            distMat = distMat + distMat';            
        else
            for i = 1 : numTestSamps
                currX = X(:,i);
                waitbar(i/numTestSamps,h)
                parfor k = 1 : numTrainSamps
                    k
                    distMat(i,k) = wde2DSlidingDist(...
                                 currX, Y(:,k),...
                                 wName, startLevel, ...
                                 stopLevel, coeffsIdx, sampleDataSupp,...
                                 slrParam,...
                        [costMatFoldPath, 'costMatrix', num2str(k), '.cm'],...
                        [assignMatFoldPath, 'assignments', num2str(k), '.txt'],...
                        execPath);
                end
            end
        end % if (isXEqualY == 1) % Distance matrix is symmetric.
end % switch distType
                             
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F);
                             
                             
                             
                             
                             
                             