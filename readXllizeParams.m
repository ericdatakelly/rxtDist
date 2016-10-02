function [textLines,rxtStructures,structureType,parameterArray...
    ] = readXllizeParams(crystallizeParamsInputFileLocation)

% get the text from the file and store it (for writing later)
% get the parameters out of the text for use now

% Get the name of the file from the calling script and then open the file
fid = fopen(crystallizeParamsInputFileLocation);

% Extract the text lines from the file so that they can be easily rewritten
% later.
textLine = fgetl(fid);
textLines = cell(0,1);
while ischar(textLine)
    textLines{end+1,1} = textLine;
    % Extract the structures into a separate array so that they can be used
    % to determine the current amount of reactant material more easily.
    if strfind(textLine,'structures')
        % Get the number of structures
        textLine = textLine(strfind(textLine,':'):length(textLine)); % find the ':' and get everything after
        numStructures = sscanf(textLine,'%*c %g'); % skip the characters and get the number
        if numStructures > 0
            % Get the structures
            rxtStructures = zeros(numStructures,8);
            textLine = fgetl(fid);
            textLine = textLine(strfind(textLine,')'):length(textLine));
            structureValues = sscanf(textLine,'%*c %g');
            structureType = structureValues(1,1);
            switch structureType
                case 1 % 1 = layer
                    rxtStructures = zeros(length(numStructures),4);
                    rxtStructures(1,:) = structureValues;
                    for n = 2:numStructures
                        textLine = fgetl(fid);
                        textLine = textLine(strfind(textLine,')'):length(textLine));
                        structureValues = sscanf(textLine,'%*c %g');
                        rxtStructures(n,:) = structureValues;
                    end
                case 2 % 2 = block
                    rxtStructures(1,:) = structureValues;
                    for n = 2:numStructures
                        textLine = fgetl(fid);
                        textLine = textLine(strfind(textLine,')'):length(textLine));
                        structureValues = sscanf(textLine,'%*c %g');
                        rxtStructures(n,:) = structureValues;
                    end
                case 3 % 3 = ellipsoid
                    rxtStructures(1,:) = structureValues;
                    for n = 2:numStructures
                        textLine = fgetl(fid);
                        textLine = textLine(strfind(textLine,')'):length(textLine));
                        structureValues = sscanf(textLine,'%*c %g');
                        rxtStructures(n,:) = structureValues;
                    end
            end
        else
            rxtStructures = zeros(1,8);
        end
    end
    textLine = fgetl(fid);
end

% Extract the values of the parameters from the text lines so that they can
% be used for calculations more easily.
textLinesSize = size(textLines);
parameterArray = zeros(textLinesSize(1,1),1);
for i = 4:length(parameterArray) % Skip headers by starting at 4
    textLine = textLines{i,1};
    textLine = textLine(strfind(textLine,':'):length(textLine)); % find the ':' and get everything after
    parameter = sscanf(textLine,'%*c %g'); % skip the characters and get the number
    if (strfind(textLines{i,1},'DelGrxn path') > 0)
        textLine = textLines{i,1};
        textLine = textLine(strfind(textLine,':'):length(textLine)); % find the ':' and get everything after
        heatingPath = sscanf(textLine,'%*c (%g,%g,%g)');
        parameterArray = [parameterArray; heatingPath]; %#ok<AGROW> % The heating path is only added once so the Matlab warning should be ignored
%         heatingPath = heatingPath';
    elseif isempty(parameter)
        parameterArray(i,1) = 0;
    else
        parameterArray(i,1) = parameter;
    end
end
fclose(fid);