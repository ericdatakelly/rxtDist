function [hetForParams...
    ] = calcLayersRandom(xllizeParamMaxX,...
    xllizeParamMaxY,xllizeParamMaxZ,heterogeneityLengthX,...
    fixedRxtConc,voxelVol,hetRandomSeed,xllizeRxtAmount,...
    xllizeParamDimensions)

disp('Random Layer Model')

rndSeed = hetRandomSeed; % Preserve the random seed

% Determine padding needed for layer thickness of single voxel
% width, odd width, or even width
if heterogeneityLengthX == 1
    heterogeneityLengthXHalf = 0;
elseif rem(heterogeneityLengthX,2) % test for odd number
    heterogeneityLengthXHalf = (heterogeneityLengthX - 1)/2;
    paddingType = 'odd';
else % other numbers must be even
    heterogeneityLengthXHalf = (heterogeneityLengthX)/2;
    paddingType = 'even';
end

% Set the amount of reactant in the starting model array to zero
rxtAmount = 0;

% Find a starting random seed
if rndSeed == 0 % If none specified, use a new random seed
    j = round(random('Uniform',1,xllizeParamMaxX));
else
    rng(rndSeed); % Tell Matlab to use the specified random seed
    j = round(random('Uniform',1,xllizeParamMaxX));
end

switch xllizeParamDimensions
    case 1 % one-dimensional model
        
        modelArrayTemp = zeros(xllizeParamMaxX);
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        while rxtAmount < xllizeRxtAmount
            modelArray = modelArrayTemp; % Use this to help stop the while loop before the max is exceeded
            % Assign ones to the center of the layer
            modelArrayTemp(j) = 1;
            thicknessX = 1; % Track the thickness of each layer for use below
            if heterogeneityLengthX > 1 % Assign ones to the rest of the layer
                switch paddingType
                    
                    case 'odd'
                        for i = 1:heterogeneityLengthXHalf
                            if j-i >= 1
                                modelArrayTemp(j-i) = 1;
                            end
                            if j+i <= xllizeParamMaxX
                                modelArrayTemp(j+i) = 1;
                            end
                        end
                        
                    case 'even'
                        for i = 1:heterogeneityLengthXHalf
                            if j-i >= 1
                                modelArrayTemp(j-i) = 1;
                                thicknessX = thicknessX + 1;
                                % If the layer is fully assigned, end the loop
                                if thicknessX == heterogeneityLengthX
                                    break
                                end
                            end
                            if j+i <= xllizeParamMaxX
                                modelArrayTemp(j+i) = 1;
                                thicknessX = thicknessX + 1;
                                % If the layer is fully assigned, end the loop
                                if thicknessX == heterogeneityLengthX
                                    break
                                end
                            end
                        end
                end
            end
            %determine the new totalRxt amount
            numVoxels = sum(modelArrayTemp);
            rxtAmount = numVoxels * voxelVol * fixedRxtConc;
            j = round(random('Uniform',1,xllizeParamMaxX));
        end
        
    case 2 % two-dimensional model
        
        modelArrayTemp = zeros(xllizeParamMaxY,xllizeParamMaxX);
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        while rxtAmount < xllizeRxtAmount
            modelArray = modelArrayTemp; % Use this to help stop the while loop before the max is exceeded
            % Assign ones to the center of the layer
            for k = 1:xllizeParamMaxY
                modelArrayTemp(k,j) = 1;
            end
            thicknessX = 1; % Track the thickness of each layer for use below
            if heterogeneityLengthX > 1 % Assign ones to the rest of the layer
                switch paddingType
                    
                    case 'odd'
                        for i = 1:heterogeneityLengthXHalf
                            if j-i >= 1
                                for k = 1:xllizeParamMaxY
                                    modelArrayTemp(k,j-i) = 1;
                                end
                            end
                            if j+i <= xllizeParamMaxX
                                for k = 1:xllizeParamMaxY
                                    modelArrayTemp(k,j+i) = 1;
                                end
                            end
                        end
                        
                    case'even'
                        for i = 1:heterogeneityLengthXHalf
                            if j-i >= 1
                                for k = 1:xllizeParamMaxY
                                    modelArrayTemp(k,j-i) = 1;
                                end
                                thicknessX = thicknessX + 1;
                                % If the layer is fully assigned, end the loop
                                if thicknessX == heterogeneityLengthX
                                    break
                                end
                            end
                            if j+i <= xllizeParamMaxX
                                for k = 1:xllizeParamMaxY
                                    modelArrayTemp(k,j+i) = 1;
                                end
                                thicknessX = thicknessX + 1;
                                % If the layer is fully assigned, end the loop
                                if thicknessX == heterogeneityLengthX
                                    break
                                end
                            end
                        end
                end
            end
            %determine the new totalRxt amount
            numVoxels = sum(sum(modelArrayTemp));
            rxtAmount = numVoxels * voxelVol * fixedRxtConc;
            j = round(random('Uniform',1,xllizeParamMaxX));
        end
        
    case 3 % three-dimensional model
        modelArrayTemp = zeros(xllizeParamMaxY,xllizeParamMaxX,...
            xllizeParamMaxZ);
        
        
        
        % Predetemine the maximum amount of reactant that can be put into
        % the model to achieve a fit.  The amount will be +/- one layer
        % thickness
        
        
        
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,2);
        [startIndexTemp,endIndexTemp] = deal(0);
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        while rxtAmount < xllizeRxtAmount
            modelArray = modelArrayTemp; % Use this to help stop the while loop before the max is exceeded
            hetForParams = [hetForParams; startIndexTemp endIndexTemp];
            
            % Assign ones to the center of the layer
            for k = 1:xllizeParamMaxY
                for m = 1:xllizeParamMaxZ
                    modelArrayTemp(k,j,m) = 1;
                end
            end
            
            thicknessX = 1; % Track the thickness of each layer for use below
            
            % Assign ones to the rest of the layer
            if heterogeneityLengthX > 1
                switch paddingType
                    
                    case 'odd'
                        for i = 1:heterogeneityLengthXHalf
                            if j-i >= 1
                                for k = 1:xllizeParamMaxY
                                    for m = 1:xllizeParamMaxZ
                                        modelArrayTemp(k,j-i,m) = 1;
                                    end
                                end
                                startIndexTemp = j-i;
                            end
                            if j+i <= xllizeParamMaxX
                                for k = 1:xllizeParamMaxY
                                    for m = 1:xllizeParamMaxZ
                                        modelArrayTemp(k,j+i,m) = 1;
                                    end
                                end
                                endIndexTemp = j+i;
                            end
                        end
                        
                    case'even'
                        for i = 1:heterogeneityLengthXHalf
                            if j-i >= 1
                                for k = 1:xllizeParamMaxY
                                    for m = 1:xllizeParamMaxZ
                                        modelArrayTemp(k,j-i,m) = 1;
                                    end
                                end
                                startIndexTemp = j-i;
                                thicknessX = thicknessX + 1;
                                % If the layer is fully assigned, end the loop
                                if thicknessX == heterogeneityLengthX
                                    break
                                end
                            end
                            if j+i <= xllizeParamMaxX
                                for k = 1:xllizeParamMaxY
                                    for m = 1:xllizeParamMaxZ
                                        modelArrayTemp(k,j+i,m) = 1;
                                    end
                                end
                                endIndexTemp = j+i;
                                thicknessX = thicknessX + 1;
                                % If the layer is fully assigned, end the loop
                                if thicknessX == heterogeneityLengthX
                                    break
                                end
                            end
                        end
                end
                
            end
            %determine the new totalRxt amount
            numVoxels = sum(sum(sum((modelArrayTemp))));
            rxtAmount = numVoxels * voxelVol * fixedRxtConc;
            
            
            %             if rxtAmount == xllizeRxtAmount
            %                 % stop making new layers
            %             elseif rxtAmount > xllizeRxtAmount
            %                 % add one layer thickness at a time to achieve a fit
            %                 % this might require a tolerance value because the
            %                 % reactant amount is discretized to +/- one layer
            %                 % if that's true, predetermine the reactant amount, given
            %                 the size of the layers, and use that predetermined amount
            %                 as the limit.  When the rxt amount exceeds the
            %                 predetermined amount, add one voxel thickness at a time
            %                 to achieve the result (do this in an if statement to
            %                 avoid running this until it is needed)
            %             else
            %             end
            
            
            j = round(random('Uniform',1,xllizeParamMaxX));
        end
end

hetForParams = hetForParams(3:end,:)-1;
hetForParams = sortrows(hetForParams);

% Check the heterogeneities for overlap
overlapIndices = 0;
for i = 2:length(hetForParams(:,1))
    if hetForParams(i,1) <= hetForParams(i-1,2)
        overlapIndices = [overlapIndices; i];
    end
end
overlapIndices = overlapIndices(2:end); % remove the zero

% Replace each pair of overlapping heterogeneities with one that spans both
for i = length(overlapIndices):-1:1
    if overlapIndices(i) == 2
        hetForParams(1,:) = [hetForParams(1,1) hetForParams(2,2)];
        hetForParams(2,:) = [];
    else
        hetForParams(overlapIndices(i)-1,:) =...
            [hetForParams(overlapIndices(i)-1,1)...
            hetForParams(overlapIndices(i),2)];
        hetForParams(overlapIndices(i),:) = [];
    end
end
