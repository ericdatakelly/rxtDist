function [hetForParams...
    ] = calcBlocksRandom(xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,heterogeneityLengthX,heterogeneityLengthY,...
    heterogeneityLengthZ,xllizeHetConc,xllizeParamDefaultRxtConc,...
    voxelVol,hetRandomSeed,xllizeRxtAmount,xllizeParamDimensions)

disp('Random Block Model')

rndSeed = hetRandomSeed; % Preserve the random seed
xllizeParamStructureType = 2;

% Use the amount of reactant in one voxel as a tolerance limit when
% deciding if the new array has the same amount of reactant as the
% previous Crystallize model
voxelRxtAmount = voxelVol * xllizeHetConc;


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
        
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX,...
            xllizeParamMaxZ); % For blocks, Crystallize requires X, Y, Z
        
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,6);
        [startIndexX,startIndexY,startIndexZ,...
            endIndexX,endIndexY,endIndexZ] = deal(0);
        
        % Set the amount of reactant in the starting model array to zero
        rxtAmount = 0;
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        while rxtAmount < xllizeRxtAmount
            hetForParams = [hetForParams;...
                startIndexX startIndexY startIndexZ...
                endIndexX endIndexY endIndexZ];
            
            % Determine the indices of the heterogeneity
            %add values to the random numbers to find the start and end
            %of each dimension
            startIndexX = x-heterogeneityLengthXHalf;
            startIndexY = 1;
            startIndexZ = 1;
            endIndexX = x+heterogeneityLengthXHalf;
            endIndexY = 1;
            endIndexZ = 1;
            
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
            modelArray(startIndexY:endIndexY,...
                startIndexX:endIndexX,...
                startIndexZ:endIndexZ) = 1;
            
            %check that the modelArray rxtAmt does not exceed the
            %maximum amount
            numVoxels = sum(sum(sum((modelArray))));
            rxtAmount = numVoxels * voxelVol * hetRxtConc;
            
            % Find new random points
            x = round(random('Uniform',1,xllizeParamMaxX));
        end
        
        % Remove leading zeros and sort the heterogeneities
        hetForParams = hetForParams(3:end,:)-1;
        hetForParams = sortrows(hetForParams);
        
        % Heterogeneities are not checked for overlap because the process
        % is too complex (I would need to divide overlapping blocks into
        % smaller individual blocks).  There's no advantage.
        
    case 2 % two-dimensional model
        
        % Determine padding needed for heterogeneity thickness
        if heterogeneityLengthX == 1
            heterogeneityLengthXHalf = 0;
        elseif rem(heterogeneityLengthX,2) % test for odd number
            heterogeneityLengthXHalf = (heterogeneityLengthX - 1)/2;
        else % other numbers must be even
            heterogeneityLengthXHalf = (heterogeneityLengthX)/2;
        end
        
        if heterogeneityLengthY == 1
            heterogeneityLengthYHalf = 0;
        elseif rem(heterogeneityLengthY,2) % test for odd number
            heterogeneityLengthYHalf = (heterogeneityLengthY - 1)/2;
        else % other numbers must be even
            heterogeneityLengthYHalf = (heterogeneityLengthY)/2;
        end
        
        % Find a starting random seed
        if rndSeed == 0 % If none specified, use a new random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
            y = round(random('Uniform',1,xllizeParamMaxX));
        else
            rng(rndSeed); % Tell Matlab to use the specified random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
            y = round(random('Uniform',1,xllizeParamMaxX));
        end
        
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX,...
            xllizeParamMaxZ); % For blocks, Crystallize requires X, Y, Z
        
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,6);
        [startIndexX,startIndexY,startIndexZ,...
            endIndexX,endIndexY,endIndexZ] = deal(0);
        
        % Set the amount of reactant in the starting model array to zero
        rxtAmount = 0;
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        while rxtAmount < xllizeRxtAmount
            hetForParams = [hetForParams;...
                startIndexX startIndexY startIndexZ...
                endIndexX endIndexY endIndexZ];
            
            % Determine the indices of the heterogeneity
            %add values to the random numbers to find the start and end
            %of each dimension
            startIndexX = x-heterogeneityLengthXHalf;
            startIndexY = y-heterogeneityLengthYHalf;
            startIndexZ = 1;
            endIndexX = x+heterogeneityLengthXHalf;
            endIndexY = y+heterogeneityLengthYHalf;
            endIndexZ = 1;
            
            %check to see if any dimension is outside of the model
            %space.  If so, make the new limit the sames as the model
            %limit
            if startIndexX < 1
                startIndexX = 1;
            end
            if startIndexY < 1
                startIndexY = 1;
            end
            if endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            end
            if endIndexY > xllizeParamMaxY
                endIndexY = xllizeParamMaxY;
            end
            
            %use the indices to write to the modelArray
            modelArray(startIndexY:endIndexY,...
                startIndexX:endIndexX,...
                startIndexZ:endIndexZ) = 1;
            
            %check that the modelArray rxtAmt does not exceed the
            %maximum amount
            numVoxels = sum(sum(sum((modelArray))));
            rxtAmount = numVoxels * voxelVol * hetRxtConc;
            
            % Find new random points
            x = round(random('Uniform',1,xllizeParamMaxX));
            y = round(random('Uniform',1,xllizeParamMaxX));
        end
        
        % Remove leading zeros and sort the heterogeneities
        hetForParams = hetForParams(3:end,:)-1;
        hetForParams = sortrows(hetForParams);
        
        % Heterogeneities are not checked for overlap because the process
        % is too complex (I would need to divide overlapping blocks into
        % smaller individual blocks).  There's no advantage.
        
    case 3 % three-dimensional model
        
        % Determine padding needed for heterogeneity thickness
        if heterogeneityLengthX == 1
            heterogeneityLengthXHalf = 0;
        elseif rem(heterogeneityLengthX,2) % test for odd number
            heterogeneityLengthXHalf = (heterogeneityLengthX - 1)/2;
        else % other numbers must be even
            heterogeneityLengthXHalf = (heterogeneityLengthX)/2;
        end
        
        if heterogeneityLengthY == 1
            heterogeneityLengthYHalf = 0;
        elseif rem(heterogeneityLengthY,2) % test for odd number
            heterogeneityLengthYHalf = (heterogeneityLengthY - 1)/2;
        else % other numbers must be even
            heterogeneityLengthYHalf = (heterogeneityLengthY)/2;
        end
        
        if heterogeneityLengthZ == 1
            heterogeneityLengthZHalf = 0;
        elseif rem(heterogeneityLengthZ,2) % test for odd number
            heterogeneityLengthZHalf = (heterogeneityLengthZ - 1)/2;
        else % other numbers must be even
            heterogeneityLengthZHalf = (heterogeneityLengthZ)/2;
        end
        
        % Find a starting random seed
        if rndSeed == 0 % If none specified, use a new random seed
            x = round(random('Uniform',1,xllizeParamMaxX));
            y = round(random('Uniform',1,xllizeParamMaxY));
            z = round(random('Uniform',1,xllizeParamMaxZ));
        else
            rng(rndSeed); % Tell Matlab to use the specified random seed
            disp('Starting seed')
            x = round(random('Uniform',1,xllizeParamMaxX));
            y = round(random('Uniform',1,xllizeParamMaxY));
            z = round(random('Uniform',1,xllizeParamMaxZ));
        end
        
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,6);
        [startIndexX,startIndexY,startIndexZ,...
            endIndexX,endIndexY,endIndexZ] = deal(0);
        
        % Set the amount of reactant in the starting model array to zero
        rxtAmount = 0;
        
        % Use this to check for a nearly full model space
        oneLayerVoxels = ones(xllizeParamMaxY,heterogeneityLengthX,...
            xllizeParamMaxZ);
        oneLayerRxtAmount = sum(sum(sum(oneLayerVoxels))) * voxelVol * xllizeHetConc;
        xllizeRxtAmountLow = xllizeRxtAmount - oneLayerRxtAmount;
        
        % Assign values to the model array and test to see if the reactant amount
        % is enough to reach the user amount or the maximum the array can
        % hold.
        % ********* This while statement is too basic.  It will eventually
        % reduce the reactants in the model significantly if it is applied
        % to the same model multiple times.
        while (rxtAmount + voxelRxtAmount) < xllizeRxtAmount
            hetForParams = [hetForParams;...
                startIndexX startIndexY startIndexZ...
                endIndexX endIndexY endIndexZ];
            
            % Determine the indices of the heterogeneity
            %add values to the random numbers to find the start and end
            %of each dimension
            startIndexX = x-heterogeneityLengthXHalf;
            startIndexY = y-heterogeneityLengthYHalf;
            startIndexZ = z-heterogeneityLengthZHalf;
            endIndexX = x+heterogeneityLengthXHalf;
            endIndexY = y+heterogeneityLengthYHalf;
            endIndexZ = z+heterogeneityLengthZHalf;
            
            %check to see if any dimension is outside of the model
            %space.  If so, make the new limit the sames as the model
            %limit
            if startIndexX < 1
                startIndexX = 1;
            end
            if startIndexY < 1
                startIndexY = 1;
            end
            if startIndexZ < 1
                startIndexZ = 1;
            end
            if endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            end
            if endIndexY > xllizeParamMaxY
                endIndexY = xllizeParamMaxY;
            end
            if endIndexZ > xllizeParamMaxZ
                endIndexZ = xllizeParamMaxZ;
            end
            
            hetForParamsTemp = [hetForParams;...
                startIndexX startIndexY startIndexZ...
                endIndexX endIndexY endIndexZ];
            %check that the modelArray rxtAmt does not exceed the
            %maximum amount
            [rxtAmount...
                ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
                xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
                xllizeParamStructureType,hetForParamsTemp,xllizeHetConc)

            % if the rxtAmt is within one block volume of the maximum,
            % add progressively smaller blocks until the maximum is reached
            while (xllizeRxtAmountLow < (rxtAmount + voxelRxtAmount)) &&...
                    ((rxtAmount + voxelRxtAmount) < xllizeRxtAmount)
                
                % Find new random numbers
                x = round(random('Uniform',1,xllizeParamMaxX));
                y = round(random('Uniform',1,xllizeParamMaxY));
                z = round(random('Uniform',1,xllizeParamMaxZ));
                
                startIndexX = x;
                startIndexY = y;
                startIndexZ = z;
                endIndexX = x + heterogeneityLengthX;
                endIndexY = y + heterogeneityLengthY;
                endIndexZ = z + heterogeneityLengthZ;
                
                if endIndexX > xllizeParamMaxX
                    endIndexX = xllizeParamMaxX;
                end
                if endIndexY > xllizeParamMaxY
                    endIndexY = xllizeParamMaxY;
                end
                if endIndexZ > xllizeParamMaxZ
                    endIndexZ = xllizeParamMaxZ;
                end
                
                hetForParamsTemp = [hetForParams;...
                    startIndexX startIndexY startIndexZ...
                    endIndexX endIndexY endIndexZ];
                
                % Check to see if enough reactant is in the new model
                [rxtAmount...
                    ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
                    xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
                    xllizeParamStructureType,hetForParamsTemp,xllizeHetConc)
                
                while (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
                    if endIndexX > startIndexX
                        endIndexX = endIndexX - 1;
                    end
                    if endIndexY > startIndexY
                        endIndexY = endIndexY - 1;
                    end
                    if endIndexY > startIndexY
                        endIndexZ = endIndexZ - 1;
                    end
                    
                    hetForParamsTemp = [hetForParams;...
                        startIndexX startIndexY startIndexZ...
                        endIndexX endIndexY endIndexZ];
                    
                    % Check to see if enough reactant is in the new model
                    [rxtAmount...
                        ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
                        xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
                        xllizeParamStructureType,hetForParamsTemp,xllizeHetConc)
                    
                    % If block is only one voxel in size, stop looping
                    if (endIndexX == startIndexX) && (endIndexY == startIndexY) && (endIndexZ == startIndexZ)
                        break
                    end
                end
                
                hetForParams = [hetForParams;...
                    startIndexX startIndexY startIndexZ...
                    endIndexX endIndexY endIndexZ];
                
            end
            if (rxtAmount + voxelRxtAmount) > xllizeRxtAmount
                %                 disp('rxtAmount >= xllizeRxtAmount #7')
                break
            end
            
            % Find new random points
            x = round(random('Uniform',1,xllizeParamMaxX));
            y = round(random('Uniform',1,xllizeParamMaxY));
            z = round(random('Uniform',1,xllizeParamMaxZ));
        end
        
        % Remove leading zeros and sort the heterogeneities
        hetForParams = hetForParams(3:end,:)-1;
        hetForParams = sortrows(hetForParams,[1 2 3]);
        
        % Report the final reactant amount for testing purposes
        [finalRxtAmount...
            ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
            xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
            xllizeParamStructureType,hetForParamsTemp,xllizeHetConc)
        
        % Heterogeneities are not checked for overlap because the process
        % is too complex (I would need to divide overlapping blocks into
        % smaller individual blocks).  There's no advantage.
end