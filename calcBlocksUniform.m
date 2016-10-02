function [hetForParams...
    ] = calcBlocksUniform(xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,gapX,gapY,gapZ,heterogeneityLengthX,...
    heterogeneityLengthY,heterogeneityLengthZ,xllizeParamDimensions)

disp('Uniform Block Model')

switch xllizeParamDimensions
    case 1 % one-dimensional model
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,6);
        
        % Loop through the matrix and add each heterogeneity
        startIndexZ = 1;
        startIndexY = 1;
        startIndexX = 1;
        endIndexY = 1;
        endIndexZ = 1;
        endIndexX = heterogeneityLengthX;
        while endIndexX < xllizeParamMaxX
            hetForParams = [hetForParams;...
                startIndexX,startIndexY,startIndexZ,...
                endIndexX,endIndexY,endIndexZ];
            % Determine the indices of the heterogeneity
            startIndexX = endIndexX + gapX + 1;
            % (The new start index is the sum of the previous end
            % index and the gap)
            endIndexX = startIndexX + heterogeneityLengthX - 1;
            %check to see if any dimension is outside of the model
            %space.  If so, make the new limit the sames as the model
            %limit
            if endIndexX > xllizeParamMaxX
                endIndexX = xllizeParamMaxX;
            end
        end
        
    case 2 % two-dimensional model
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,6);
        
        % Loop through the matrix and add each heterogeneity
        startIndexZ = 1;
        startIndexY = 1;
        endIndexZ = 1;
        endIndexY = heterogeneityLengthY;
        while endIndexY < xllizeParamMaxY
            startIndexX = 1;
            endIndexX = heterogeneityLengthX;
            while endIndexX < xllizeParamMaxX
                hetForParams = [hetForParams;...
                    startIndexX,startIndexY,startIndexZ,...
                    endIndexX,endIndexY,endIndexZ];
                % Determine the indices of the heterogeneity
                startIndexX = endIndexX + gapX + 1;
                % (The new start index is the sum of the previous end
                % index and the gap)
                endIndexX = startIndexX + heterogeneityLengthX - 1;
                %check to see if any dimension is outside of the model
                %space.  If so, make the new limit the sames as the model
                %limit
                if endIndexX > xllizeParamMaxX
                    endIndexX = xllizeParamMaxX;
                end
            end
            startIndexY = endIndexY + gapY + 1;
            endIndexY = startIndexY + heterogeneityLengthY - 1;
            if endIndexY > xllizeParamMaxY
                endIndexY = xllizeParamMaxY;
            end
        end
        
    case 3 % three-dimensional model
        % Record the heterogeneities for writing to the params file
        hetForParams = zeros(1,6);
        
        % Loop through the matrix and add each heterogeneity
        startIndexZ = 1;
        endIndexZ = heterogeneityLengthZ;
        while endIndexZ < xllizeParamMaxZ
            startIndexY = 1;
            endIndexY = heterogeneityLengthY;
            while endIndexY < xllizeParamMaxY
                startIndexX = 1;
                endIndexX = heterogeneityLengthX;
                while endIndexX < xllizeParamMaxX
                    hetForParams = [hetForParams;...
                        startIndexX,startIndexY,startIndexZ,...
                        endIndexX,endIndexY,endIndexZ];
                    % Determine the indices of the heterogeneity
                    startIndexX = endIndexX + gapX + 1;
                    % (The new start index is the sum of the previous end
                    % index and the gap)
                    endIndexX = startIndexX + heterogeneityLengthX - 1;
                    %check to see if any dimension is outside of the model
                    %space.  If so, make the new limit the sames as the model
                    %limit
                    if endIndexX > xllizeParamMaxX
                        endIndexX = xllizeParamMaxX;
                    end
                end
                startIndexY = endIndexY + gapY + 1;
                endIndexY = startIndexY + heterogeneityLengthY - 1;
                if endIndexY > xllizeParamMaxY
                    endIndexY = xllizeParamMaxY;
                end
            end
            startIndexZ = endIndexZ + gapZ + 1;
            endIndexZ = startIndexZ + heterogeneityLengthZ - 1;
            if endIndexZ > xllizeParamMaxZ
                endIndexZ = xllizeParamMaxZ;
            end
        end
        
        % Remove leading zeros and sort the heterogeneities
        hetForParams = hetForParams(2:end,:)-1;
        hetForParams = sortrows(hetForParams);
end
