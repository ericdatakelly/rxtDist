function [hetForParams...
    ] = calcLayersRandom(xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,heterogeneityLengthX,xllizeHetConc,...
    xllizeParamDefaultRxtConc,voxelVol,hetRandomSeed,...
    xllizeRxtAmount,xllizeParamDimensions)

disp('Random Layer Model')

rndSeed = hetRandomSeed; % Preserve the random seed
xllizeParamStructureType = 1;

switch xllizeParamDimensions
    case 1 % one-dimensional model
        
        % Determine padding needed for heterogeneity thickness
        if heterogeneityLengthX == 1
            heterogeneityLengthXHalf = 0;
        elseif rem(heterogeneityLengthX,2) % test for odd number
            heterogeneityLengthXHalf = (heterogeneityLengthX - 1)/2;
        else % other numbers must be even
            heterogeneityLengthXHalf = (heterogeneityLengthX)/2;
        end
        
        % Find a starting random seed
        if rndSeed == 0 % If none specified, use a new random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
        else
            rng(rndSeed); % Tell Matlab to use the specified random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
        end
        
        modelArray = zeros(xllizeParamMaxX);
        
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,2);
        [startIndexX,endIndexX] = deal(0);
        
        % Set the amount of reactant in the starting model array to zero
        rxtAmount = 0;
        
        % Use the amount of reactant in one voxel as a tolerance limit when
        % deciding if the new array has the same amount of reactant as the
        % previous Crystallize model
        voxelRxtAmount = voxelVol * xllizeHetConc;
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        % Determine the amount of reactant in one layer and then check to
        % see if rxt amount of model is within one layer of the max.
        oneLayerVoxels = ones(heterogeneityLengthX);
        oneLayerRxtAmount = sum(oneLayerVoxels) * voxelVol * xllizeHetConc;
        xllizeRxtAmountLow = xllizeRxtAmount - oneLayerRxtAmount;
        while (rxtAmount + voxelRxtAmount) < xllizeRxtAmount
            hetForParams = [hetForParams; startIndexX endIndexX];
            
            % Determine the indices of the heterogeneity
            %add values to the random numbers to find the start and end
            %of each dimension
            startIndexX = x-heterogeneityLengthXHalf;
            endIndexX = x+heterogeneityLengthXHalf;
            
            %check to see if any dimension is outside of the model
            %space.  If so, make the new limit the sames as the model
            %limit
            if startIndexX < 1
                startIndexX = 1;
            end
            if endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            end
            
            %use the indices to write to the modelArray
%             disp('Write to array.  Still less than max.')
            modelArray(startIndexX:endIndexX) = 1;
            
            %check that the modelArray rxtAmt does not exceed the
            %maximum amount
            numVoxels = sum(modelArray);
            rxtAmount = numVoxels * voxelVol * xllizeHetConc;
            
            % if the rxtAmt is within one layer thickness of the maximum,
            % add progressively larger layers until the maximum is reached
            while (xllizeRxtAmountLow < (rxtAmount + voxelRxtAmount)) &&...
                    ((rxtAmount + voxelRxtAmount) < xllizeRxtAmount)
                % Find a new random number
                x = round(random('Uniform',1,xllizeParamMaxX));
                startIndexX = x;
                endIndexX = x;
                % Write to the array at the new coordinate
%                 fprintf('Amount between low and max.  Try %i\n',j)
                modelArray(x) = 1;
                % Check to see if enough reactant is in the new model
                numVoxels = sum(modelArray);
                rxtAmount = numVoxels * voxelVol * xllizeHetConc;
                if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                     disp('rxtAmount >= xllizeRxtAmount #1')
%                     hetForParams = [hetForParams; startIndexX endIndexX]
                    break
                end
                % Write single-voxel widths to the layer until the reactant
                % amount has been reached
                for i = 1:heterogeneityLengthXHalf
                    if x-i >= 1
%                         fprintf('Try startIndex %i\n',j-i)
                        modelArray(x-i) = 1;
                        % Check to see if enough reactant is in the new model
                        numVoxels = sum(modelArray);
                        rxtAmount = numVoxels * voxelVol * xllizeHetConc;
                        if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                             disp('rxtAmount >= xllizeRxtAmount #3')
                            break
                        end
                        startIndexX = x-i;
                    end
                    if x+i <= xllizeParamMaxX
%                         fprintf('Try endIndex %i\n',j+i)
                        modelArray(x+i) = 1;
                        % Check to see if enough reactant is in the new model
                        numVoxels = sum(modelArray);
                        rxtAmount = numVoxels * voxelVol * xllizeHetConc;
                        if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                             disp('rxtAmount >= xllizeRxtAmount #5')
                            break
                        end
                        endIndexX = x+i;
                    end
                end
%                 fprintf('Final indices for one het: %i\t%i\n',...
%                     startIndexX,endIndexX)
                hetForParams = [hetForParams; startIndexX endIndexX];
                if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                     disp('rxtAmount >= xllizeRxtAmount #6')
                    break
                end
            end
            if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                 disp('rxtAmount >= xllizeRxtAmount #7')
                break
            end
            % Find new random points
            x = round(random('Uniform',1,xllizeParamMaxX));
        end
        
    case 2 % two-dimensional model
        
        % Determine padding needed for heterogeneity thickness
        if heterogeneityLengthX == 1
            heterogeneityLengthXHalf = 0;
        elseif rem(heterogeneityLengthX,2) % test for odd number
            heterogeneityLengthXHalf = (heterogeneityLengthX - 1)/2;
        else % other numbers must be even
            heterogeneityLengthXHalf = (heterogeneityLengthX)/2;
        end
        
        % Find a starting random seed
        if rndSeed == 0 % If none specified, use a new random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
        else
            rng(rndSeed); % Tell Matlab to use the specified random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
        end
        
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX);
        
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,2);
        [startIndexX,endIndexX] = deal(0);
        
        % Set the amount of reactant in the starting model array to zero
        rxtAmount = 0;
        
        % Use the amount of reactant in one voxel as a tolerance limit when
        % deciding if the new array has the same amount of reactant as the
        % previous Crystallize model
        voxelRxtAmount = voxelVol * xllizeHetConc;
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        % Determine the amount of reactant in one layer and then check to
        % see if rxt amount of model is within one layer of the max.
        oneLayerVoxels = ones(xllizeParamMaxY,heterogeneityLengthX);
        oneLayerRxtAmount = sum(sum(oneLayerVoxels)) * voxelVol * xllizeHetConc;
        xllizeRxtAmountLow = xllizeRxtAmount - oneLayerRxtAmount;
        while (rxtAmount + voxelRxtAmount) < xllizeRxtAmount
            hetForParams = [hetForParams; startIndexX endIndexX];
            
            % Determine the indices of the heterogeneity
            %add values to the random numbers to find the start and end
            %of each dimension
            startIndexX = x-heterogeneityLengthXHalf;
            endIndexX = x+heterogeneityLengthXHalf;
            
            %check to see if any dimension is outside of the model
            %space.  If so, make the new limit the sames as the model
            %limit
            if startIndexX < 1
                startIndexX = 1;
            end
            if endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            end
            
            %use the indices to write to the modelArray
%             disp('Write to array.  Still less than max.')
            modelArray(1:xllizeParamMaxY,startIndexX:endIndexX) = 1;
            
            %check that the modelArray rxtAmt does not exceed the
            %maximum amount
            numVoxels = sum(sum(sum(modelArray)));
            rxtAmount = numVoxels * voxelVol * xllizeHetConc;
            
            % if the rxtAmt is within one layer thickness of the maximum,
            % add progressively larger layers until the maximum is reached
            while (xllizeRxtAmountLow < (rxtAmount + voxelRxtAmount)) &&...
                    ((rxtAmount + voxelRxtAmount) < xllizeRxtAmount)
                % Find a new random number
                x = round(random('Uniform',1,xllizeParamMaxX));
                startIndexX = x;
                endIndexX = x;
                % Write to the array at the new coordinate
%                 fprintf('Amount between low and max.  Try %i\n',j)
                modelArray(1:xllizeParamMaxY,x) = 1;
                % Check to see if enough reactant is in the new model
                numVoxels = sum(sum(modelArray));
                rxtAmount = numVoxels * voxelVol * xllizeHetConc;
                if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                     disp('rxtAmount >= xllizeRxtAmount #1')
%                     hetForParams = [hetForParams; startIndexX endIndexX]
                    break
                end
                % Write single-voxel widths to the layer until the reactant
                % amount has been reached
                for i = 1:heterogeneityLengthXHalf
                    if x-i >= 1
%                         fprintf('Try startIndex %i\n',j-i)
                        modelArray(1:xllizeParamMaxY,x-i) = 1;
                        % Check to see if enough reactant is in the new model
                        numVoxels = sum(sum(modelArray));
                        rxtAmount = numVoxels * voxelVol * xllizeHetConc;
                        if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                             disp('rxtAmount >= xllizeRxtAmount #3')
                            break
                        end
                        startIndexX = x-i;
                    end
                    if x+i <= xllizeParamMaxX
%                         fprintf('Try endIndex %i\n',j+i)
                        modelArray(1:xllizeParamMaxY,x+i) = 1;
                        % Check to see if enough reactant is in the new model
                        numVoxels = sum(sum(modelArray));
                        rxtAmount = numVoxels * voxelVol * xllizeHetConc;
                        if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                             disp('rxtAmount >= xllizeRxtAmount #5')
                            break
                        end
                        endIndexX = x+i;
                    end
                end
%                 fprintf('Final indices for one het: %i\t%i\n',...
%                     startIndexX,endIndexX)
                hetForParams = [hetForParams; startIndexX endIndexX];
                if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                     disp('rxtAmount >= xllizeRxtAmount #6')
                    break
                end
            end
            if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                 disp('rxtAmount >= xllizeRxtAmount #7')
                break
            end
            % Find new random points
            x = round(random('Uniform',1,xllizeParamMaxX));
        end
        
    case 3 % three-dimensional model
        
        % Determine padding needed for heterogeneity thickness
        if heterogeneityLengthX == 1
            heterogeneityLengthXHalf = 0;
        elseif rem(heterogeneityLengthX,2) % test for odd number
            heterogeneityLengthXHalf = (heterogeneityLengthX - 1)/2;
        else % other numbers must be even
            heterogeneityLengthXHalf = (heterogeneityLengthX)/2;
        end
        
        % Find a starting random seed
        if rndSeed == 0 % If none specified, use a new random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
        else
            rng(rndSeed); % Tell Matlab to use the specified random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
        end
        
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX,...
            xllizeParamMaxZ);
        
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,2);
        [startIndexX,endIndexX] = deal(0);
        
        % Determine the initial reactant amount
        [rxtAmount...
            ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
            xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
            xllizeParamStructureType,hetForParams,xllizeHetConc);
        
        % Use the amount of reactant in one voxel as a tolerance limit when
        % deciding if the new array has the same amount of reactant as the
        % previous Crystallize model
        voxelRxtAmount = voxelVol * xllizeHetConc;
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        % Determine the amount of reactant in one layer and then check to
        % see if rxt amount of model is within one layer of the max.
        oneLayerVoxels = ones(xllizeParamMaxY,heterogeneityLengthX,...
            xllizeParamMaxZ);
        oneLayerRxtAmount = sum(sum(sum(oneLayerVoxels))) * voxelVol * xllizeHetConc;
        xllizeRxtAmountLow = xllizeRxtAmount - oneLayerRxtAmount;
        while (rxtAmount + voxelRxtAmount) < xllizeRxtAmount
            hetForParams = [hetForParams; startIndexX endIndexX];
            
            % Determine the indices of the heterogeneity
            %add values to the random numbers to find the start and end
            %of each dimension
            startIndexX = x-heterogeneityLengthXHalf;
            endIndexX = x+heterogeneityLengthXHalf;
            
            %check to see if any dimension is outside of the model
            %space.  If so, make the new limit the sames as the model
            %limit
            if startIndexX < 1
                startIndexX = 1;
            end
            if endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            end
            
            %use the indices to write to the modelArray
%             disp('Write to array.  Still less than max.')
            modelArray(1:xllizeParamMaxY,startIndexX:endIndexX,...
                1:xllizeParamMaxZ) = 1;
            
            %check that the modelArray rxtAmt does not exceed the
            %maximum amount
            [rxtAmount...
                ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
                xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
                xllizeParamStructureType,hetForParams,xllizeHetConc);
            
            % if the rxtAmt is within one layer thickness of the maximum,
            % add progressively larger layers until the maximum is reached
            while (xllizeRxtAmountLow < (rxtAmount + voxelRxtAmount)) &&...
                    ((rxtAmount + voxelRxtAmount) < xllizeRxtAmount)
                % Find a new random number
                x = round(random('Uniform',1,xllizeParamMaxX));
                startIndexX = x;
                endIndexX = x;
                % Write to the array at the new coordinate
%                 fprintf('Amount between low and max.  Try %i\n',j)
                modelArray(1:xllizeParamMaxY,x,1:xllizeParamMaxZ) = 1;
                % Check to see if enough reactant is in the new model
                [rxtAmount...
                    ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
                    xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
                    xllizeParamStructureType,hetForParams,xllizeHetConc);
                if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                     disp('rxtAmount >= xllizeRxtAmount #1')
%                     hetForParams = [hetForParams; startIndexX endIndexX]
                    break
                end
                % Write single-voxel widths to the layer until the reactant
                % amount has been reached
                for i = 1:heterogeneityLengthXHalf
                    if x-i >= 1
%                         fprintf('Try startIndex %i\n',j-i)
                        modelArray(1:xllizeParamMaxY,x-i,1:xllizeParamMaxZ) = 1;
                        % Check to see if enough reactant is in the new model
                        [rxtAmount...
                            ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
                            xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
                            xllizeParamStructureType,hetForParams,xllizeHetConc);
                        if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                             disp('rxtAmount >= xllizeRxtAmount #3')
                            break
                        end
                        startIndexX = x-i;
                    end
                    if x+i <= xllizeParamMaxX
%                         fprintf('Try endIndex %i\n',j+i)
                        modelArray(1:xllizeParamMaxY,x+i,1:xllizeParamMaxZ) = 1;
                        % Check to see if enough reactant is in the new model
                        [rxtAmount...
                            ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
                            xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
                            xllizeParamStructureType,hetForParams,xllizeHetConc);
                        if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                             disp('rxtAmount >= xllizeRxtAmount #5')
                            break
                        end
                        endIndexX = x+i;
                    end
                end
%                 fprintf('Final indices for one het: %i\t%i\n',...
%                     startIndexX,endIndexX)
                hetForParams = [hetForParams; startIndexX endIndexX];
                if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                     disp('rxtAmount >= xllizeRxtAmount #6')
                    break
                end
            end
            if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
%                 disp('rxtAmount >= xllizeRxtAmount #7')
                break
            end
            % Find new random points
            x = round(random('Uniform',1,xllizeParamMaxX));
        end
end

% Remove leading zeros and sort the heterogeneities
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
