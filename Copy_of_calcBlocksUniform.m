function [modelArray...
    ] = calcBlocksUniform(xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,fixedGapX,heterogeneityLengthX,...
    xllizeParamDimensions)

% rxtDistLayersUniform.m
% syntax: rxtDistLayersUniform()
% This script creates uniformly distributed blocks of reactant 
% concentration for use in Crystallize3D (Ketcham and Carlson, 2012)
% More description...

% Description (needs editing): User chooses the model name, reactant block dimensions, and reactant
% concentration.  The script distributes blocks with the same dimensions.
% If Crystallize is given a default reactant concentration, this script
% will only write over the new layers and leave the gaps between the layers
% at the default concentration.  Heterogeneities must be sized in voxel 
% units (no fractional voxels sizes)

disp('Uniform Block Model')

switch xllizeParamDimensions
    
    case 1 % one-dimensional model
        % Create starting array and assign indices to the first heterogeneity
        modelArray = zeros(1,xllizeParamMaxX);
        startIndexX = 1;
        endIndexX = heterogeneityLengthX;
        % Overwrite the array at the heterogeneity coordinates (later give the
        % array concentration values)
        while endIndexX < xllizeParamMaxX
            % Write the first heterogeneity to the array
            for j = startIndexX:endIndexX
                modelArray(j) = 1;
            end
            % Increment the starting indices to write the next heterogeneity
            startIndexX = endIndexX + fixedLayerGapX + 1;
            endIndexX = endIndexX + fixedLayerGapX + heterogeneityLengthX;
            if startIndexX > xllizeParamMaxX
                break % Stop adding heterogeneities
            elseif startIndexX == xllizeParamMaxX
                endIndexX = startIndexX; % The last heterogeneity will be only one voxel thickness
            elseif endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            else % Keep the new layerEndIndex value and write the heterogeneity
            end
        end
    case 2 % two-dimensional model
        % Create starting array and assign indices to the first heterogeneity
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX);
        [startIndexX,startIndexY] = deal(1);
        endIndexX = heterogeneityLengthX;
        endIndexY = heterogeneityLengthY;
        % Overwrite the array at the heterogeneity coordinates (later give the
        % array concentration values)
        while endIndexX < xllizeParamMaxX
            % Write the first heterogeneity to the array
            for j = startIndexX:endIndexX
                while endIndexY < xllizeParamMaxY
                    for k = startIndexY:endIndexY
                        modelArray(k,j) = 1;
                    end
                end
                % Increment the starting indices to write the next heterogeneity
                startIndexY = endIndexY + fixedLayerGapY + 1;
                endIndexY = endIndexY + fixedLayerGapY + heterogeneityLengthY;
                if startIndexY > xllizeParamMaxY
                    break % Stop adding heterogeneities
                elseif startIndexY == xllizeParamMaxY
                    endIndexY = startIndexY; % The last heterogeneity will be only one voxel thickness
                elseif endIndexY > xllizeParamMaxY
                    endIndexY = xllizeParamMaxY;
                else % Keep the new endIndex value and write the heterogeneity
                end
            end
            
            % Increment the starting indices to write the next heterogeneity
            startIndexX = endIndexX + fixedLayerGapX + 1;
            endIndexX = endIndexX + fixedLayerGapX + heterogeneityLengthX;
            if startIndexX > xllizeParamMaxX
                break % Stop adding heterogeneities
            elseif startIndexX == xllizeParamMaxX
                endIndexX = startIndexX; % The last heterogeneity will be only one voxel thickness
            elseif endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            else % Keep the new layerEndIndex value and write the heterogeneity
            end
        end
    case 3 % three-dimensional model
        % Create starting array and assign indices to the first heterogeneity
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX,xllizeParamMaxZ);
        [startIndexX,startIndexY,startIndexZ] = deal(1);
        endIndexX = heterogeneityLengthX;
        endIndexY = heterogeneityLengthY;
        endIndexZ = heterogeneityLengthZ;
        % Overwrite the array at the heterogeneity coordinates (later give the
        % array concentration values)
        while endIndexX < xllizeParamMaxX
            % Write the first heterogeneity to the array
            for j = startIndexX:endIndexX
                while endIndexY < xllizeParamMaxY
                    for k = startIndexY:endIndexY
                        while endInexZ < xllizeParamMaxZ
                            for m = startIndexZ:endIndexZ
                                modelArray(k,j,m) = 1;
                            end
                            % Increment the starting indices to write the next heterogeneity
                            startIndexZ = endIndexZ + fixedLayerGapZ + 1;
                            endIndexZ = endIndexZ + fixedLayerGapZ + heterogeneityLengthZ;
                            if startIndexZ > xllizeParamMaxZ
                                break % Stop adding heterogeneities
                            elseif startIndexZ == xllizeParamMaxZ
                                endIndexZ = startIndexZ; % The last heterogeneity will be only one voxel thickness
                            elseif endIndexZ > xllizeParamMaxZ
                                endIndexZ = xllizeParamMaxZ;
                            else % Keep the new layerEndIndex value and write the heterogeneity
                            end
                        end
                    end
                    % Increment the starting indices to write the next heterogeneity
                    startIndexY = endIndexY + fixedLayerGapY + 1;
                    endIndexY = endIndexY + fixedLayerGapY + heterogeneityLengthY;
                    if startIndexY > xllizeParamMaxY
                        break % Stop adding heterogeneities
                    elseif startIndexY == xllizeParamMaxY
                        endIndexY = startIndexY; % The last heterogeneity will be only one voxel thickness
                    elseif endIndexY > xllizeParamMaxY
                        endIndexY = xllizeParamMaxY;
                    else % Keep the new endIndex value and write the heterogeneity
                    end
                end
            end
            % Increment the starting indices to write the next heterogeneity
            startIndexX = endIndexX + fixedLayerGapX + 1;
            endIndexX = endIndexX + fixedLayerGapX + heterogeneityLengthX;
            if startIndexX > xllizeParamMaxX
                break % Stop adding heterogeneities
            elseif startIndexX == xllizeParamMaxX
                endIndexX = startIndexX; % The last heterogeneity will be only one voxel thickness
            elseif endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            else % Keep the new layerEndIndex value and write the heterogeneity
            end
        end
end