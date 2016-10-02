function [xllizeRxtAmount,xllizeHetConc...
    ] = calcXllizeRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
    xllizeParamStructureType,xllizeParamNumStructures,rxtStructures)

switch xllizeParamDimensions
    case 1 % one-dimensional model
        disp('Calculate 1D array - needs code!')
    case 2 % two-dimensional model
        disp('Calculate 2D array')
        if xllizeParamNumStructures == 0 % checks for homogeneous model
            % use number of dimensions, max voxels, and default reactant amount to
            % calculate grams of CAR
            xllizeRxtAmount = voxelVol * xllizeParamMaxX * xllizeParamMaxY *...
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
                    for i = 1:xllizeParamNumStructures
                        startIndexX = rxtStructures(i,2)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        endIndexX = rxtStructures(i,3)+1;
                        xllizeHetConc = rxtStructures(i,4);
                        xllizeArray(startIndexX:endIndexX,...
                            1:xllizeParamMaxY) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    xllizeRxtAmount = sum(sum(sum(xllizeArray * voxelVol)));
                    
                case 2 % 2 = block
                    % first make an array with the default amount
                    xllizeArray = ones(xllizeParamMaxX,xllizeParamMaxY);
                    xllizeArray = xllizeArray * xllizeParamDefaultRxtConc;
                    % then change array to reflect the new heteregeneities
                    for i = 1:xllizeParamNumStructures
                        startIndexX = rxtStructures(i,2)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        startIndexY = rxtStructures(i,3)+1;
                        endIndexX = rxtStructures(i,3)+1;
                        endIndexY = rxtStructures(i,4)+1;
                        xllizeHetConc = rxtStructures(i,8);
                        xllizeArray(startIndexX:endIndexX,...
                            startIndexY:endIndexY) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    xllizeRxtAmount = sum(sum(xllizeArray * voxelVol));
                    
                case 3 % 3 = ellipsoid
                    % ------------> *** READ THIS *** <-----------------------------------------------------------
                    % ask Rich how he assigns the ellipsiod coordinates so
                    % that I can reproduce the Crystallize3D code and keep
                    % the params files consistent.
            end
        end
        
    case 3 % three-dimensional model
        disp('Calculate 3D array')
        if xllizeParamNumStructures == 0 % checks for homogeneous model
            % use number of dimensions, max voxels, and default reactant amount to
            % calculate grams of CAR
            xllizeRxtAmount = voxelVol * xllizeParamMaxX * xllizeParamMaxY *...
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
                    for i = 1:xllizeParamNumStructures
                        startIndexX = rxtStructures(i,2)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        endIndexX = rxtStructures(i,3)+1;
                        xllizeHetConc = rxtStructures(i,4);
                        xllizeArray(startIndexX:endIndexX,...
                            1:xllizeParamMaxY,...
                            1:xllizeParamMaxZ) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    xllizeRxtAmount = sum(sum(sum(xllizeArray * voxelVol)));
                    
                case 2 % 2 = block
                    % first make an array with the default amount
                    xllizeArray = ones(xllizeParamMaxX,xllizeParamMaxY,...
                        xllizeParamMaxZ);
                    xllizeArray = xllizeArray * xllizeParamDefaultRxtConc;
                    % then change array to reflect the new heteregeneities
                    for i = 1:xllizeParamNumStructures
                        startIndexX = rxtStructures(i,2)+1; % Crystallize indices start at zero so these indices need to be shifted by one
                        startIndexY = rxtStructures(i,3)+1;
                        startIndexZ = rxtStructures(i,4)+1;
                        endIndexX = rxtStructures(i,5)+1;
                        endIndexY = rxtStructures(i,6)+1;
                        endIndexZ = rxtStructures(i,7)+1;
                        xllizeHetConc = rxtStructures(i,8);
                        xllizeArray(startIndexX:endIndexX,...
                            startIndexY:endIndexY,...
                            startIndexZ:endIndexZ) = xllizeHetConc;
                    end
                    % then sum the grams of reactant
                    xllizeRxtAmount = sum(sum(sum(xllizeArray * voxelVol)));
                    
                case 3 % 3 = ellipsoid
                    % ------------> *** READ THIS *** <-----------------------------------------------------------
                    % ask Rich how he assigns the ellipsiod coordinates so
                    % that I can reproduce the Crystallize3D code and keep
                    % the params files consistent.
            end
        end
end