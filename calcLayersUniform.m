function [hetForParams...
    ] = calcLayersUniform(xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,gapX,heterogeneityLengthX,...
    xllizeParamDimensions)

% rxtDistLayersUniform.m
% syntax: rxtDistLayersUniform()
% This script creates uniformly distributed layers of reactant 
% concentration for use in Crystallize3D (Ketcham and Carlson, 2012)
% More description...

% User chooses the model name, reactant layer thickness, and reactant
% concentration.  The script distributes layers with the same thickness.
% If Crystallize is given a default reactant concentration, this script
% will only write over the new layers and leave the gaps between the layers
% at the default concentration.  Heterogeneities must be sized in voxel 
% units (no half voxels, etc.)

disp('Uniform Layer Model')

switch xllizeParamDimensions
    
    case 1 % one-dimensional model
        modelArray = zeros(xllizeParamMaxX);
        startIndex = 1;
        endIndex = heterogeneityLengthX;
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,2);
        while endIndex < xllizeParamMaxX
            modelArray(startIndex:endIndex) = 1;
            hetForParams = [hetForParams;startIndex endIndex];
            startIndex = endIndex + gapX + 1;
            endIndex = startIndex + heterogeneityLengthX - 1;
            if startIndex > xllizeParamMaxX
                break % Stop adding layers
            elseif startIndex == xllizeParamMaxX
                endIndex = startIndex; % The last layer will be only one voxel thick
            elseif endIndex > xllizeParamMaxX
                endIndex = xllizeParamMaxX;
            else % Keep the new endIndex value and write the layer
            end
        end
        
    case 2 % two-dimensional model
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX);
        startIndex = 1;
        endIndex = heterogeneityLengthX;
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,2);
        while endIndex < xllizeParamMaxX
            modelArray(1:xllizeParamMaxY,startIndex:endIndex) = 1;
            hetForParams = [hetForParams;startIndex endIndex];
            startIndex = endIndex + gapX + 1;
            endIndex = startIndex + heterogeneityLengthX - 1;
            if startIndex > xllizeParamMaxX
                break % Stop adding layers
            elseif startIndex == xllizeParamMaxX
                endIndex = startIndex; % The last layer will be only one voxel thick
            elseif endIndex > xllizeParamMaxX
                endIndex = xllizeParamMaxX;
            else % Keep the new endIndex value and write the layer
            end
        end
        
    case 3 % three-dimensional model
        modelArray = zeros(xllizeParamMaxY,xllizeParamMaxX,xllizeParamMaxZ);
        startIndex = 1;
        endIndex = heterogeneityLengthX;
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,2);
        while endIndex < xllizeParamMaxX
            modelArray(1:xllizeParamMaxY,startIndex:endIndex,1:xllizeParamMaxZ) = 1;
            hetForParams = [hetForParams;startIndex endIndex];
            startIndex = endIndex + gapX + 1;
            endIndex = startIndex + heterogeneityLengthX - 1;
            if startIndex > xllizeParamMaxX
                break % Stop adding layers
            elseif startIndex == xllizeParamMaxX
                endIndex = startIndex; % The last layer will be only one voxel thick
            elseif endIndex > xllizeParamMaxX
                endIndex = xllizeParamMaxX;
            else % Keep the new endIndex value and write the layer
            end
        end
end
hetForParams = hetForParams(2:end,:)-1;