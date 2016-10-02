function [newRxtAmount...
    ] = calcNewRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
    xllizeParamStructureType,hetForParams,xllizeHetConc)

% Remove leading zeros in array
hetForParams = hetForParams(find(hetForParams(:,1),1):end,:);

newNumStructures = length(hetForParams(:,1));

switch xllizeParamDimensions
    case 1 % one-dimensional model
        disp('Calculate 1D array - needs code!')
    case 2 % two-dimensional model
        disp('Calculate 2D array')
        if newNumStructures == 0 % checks for homogeneous model
            % use number of dimensions, max voxels, and default reactant amount to
            % calculate grams of CAR
            newRxtAmount = voxelVol * xllizeParamMaxX * xllizeParamMaxY *...
                xllizeParamDefaultRxtConc;
        else
            % calculate the grams of CAR by summing the heterogeneities
            % use the format codes (layer, block, ...) to switch between methods
            % here
            switch xllizeParamStructureType
                case 1 % 1 = layer
                    % first make an array with the default amount
                    xllizeArray = ones(xllizeParamMaxX,xllizeParamMaxY);
                    xllizeArray = xllizeArray * xllizeParamDefaultRxtConc;
                    % then change array to reflect the new heteregeneities
                    for i = 1:newNumStructures
                        startIndexX = hetForParams(i,1)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        endIndexX = hetForParams(i,2)+1;
                        xllizeArray(startIndexX:endIndexX,...
                            1:xllizeParamMaxY) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    newRxtAmount = sum(sum(sum(xllizeArray * voxelVol)));
                    
                case 2 % 2 = block
                    % first make an array with the default amount
                    xllizeArray = ones(xllizeParamMaxX,xllizeParamMaxY);
                    xllizeArray = xllizeArray * xllizeParamDefaultRxtConc;
                    % then change array to reflect the new heteregeneities
                    for i = 1:newNumStructures
                        startIndexX = hetForParams(i,1)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        startIndexY = hetForParams(i,2)+1;
                        endIndexX = hetForParams(i,4)+1;
                        endIndexY = hetForParams(i,5)+1;
                        xllizeArray(startIndexX:endIndexX,...
                            startIndexY:endIndexY) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    newRxtAmount = sum(sum(xllizeArray * voxelVol));
                    
                case 3 % 3 = ellipsoid
                    % ------------> *** READ THIS *** <-----------------------------------------------------------
                    % ask Rich how he assigns the ellipsiod coordinates so
                    % that I can reproduce the Crystallize3D code and keep
                    % the params files consistent.
            end
        end
        
    case 3 % three-dimensional model
%         disp('Check 3D array')
        if newNumStructures == 0 % checks for homogeneous model
            % use number of dimensions, max voxels, and default reactant amount to
            % calculate grams of CAR
            newRxtAmount = voxelVol * xllizeParamMaxX * xllizeParamMaxY *...
                xllizeParamMaxZ * xllizeParamDefaultRxtConc;
        else
            % calculate the grams of CAR by summing the heterogeneities
            % use the format codes (layer, block, ...) to switch between methods
            % here
            switch xllizeParamStructureType
                
                case 1 % 1 = layer
                    % first make an array with the default concentration
                    xllizeArray = ones(xllizeParamMaxX,xllizeParamMaxY,...
                        xllizeParamMaxZ);
                    xllizeArray = xllizeArray * xllizeParamDefaultRxtConc;
                    % then overwrite the array with the new heteregeneities
                    for i = 1:newNumStructures
                        startIndexX = hetForParams(i,1)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        endIndexX = hetForParams(i,2)+1;
                        xllizeArray(1:xllizeParamMaxY,...
                            startIndexX:endIndexX,...
                            1:xllizeParamMaxZ) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    newRxtAmount = sum(sum(sum(xllizeArray * voxelVol)));
                    
                case 2 % 2 = block
                    % first make an array with the default amount
                    xllizeArray = ones(xllizeParamMaxY,xllizeParamMaxX,...
                        xllizeParamMaxZ);
                    xllizeArray = xllizeArray * xllizeParamDefaultRxtConc;
                    % then change array to reflect the new heteregeneities
                    for i = 1:newNumStructures
                        startIndexX = hetForParams(i,1)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        startIndexY = hetForParams(i,2)+1;
                        startIndexZ = hetForParams(i,3)+1;
                        endIndexX = hetForParams(i,4)+1;
                        endIndexY = hetForParams(i,5)+1;
                        endIndexZ = hetForParams(i,6)+1;
                        xllizeArray(startIndexY:endIndexY,...
                            startIndexX:endIndexX,...
                            startIndexZ:endIndexZ) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    newRxtAmount = sum(sum(sum(xllizeArray * voxelVol)));
                    
                case 3 % 3 = ellipsoid
                    % ------------> *** READ THIS *** <-----------------------------------------------------------
                    % ask Rich how he assigns the ellipsiod coordinates so
                    % that I can reproduce the Crystallize3D code and keep
                    % the params files consistent.
            end
        end
end