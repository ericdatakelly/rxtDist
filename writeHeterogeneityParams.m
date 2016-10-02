function writeHeterogeneityParams(heterogeneityParamsOutputFileLocation,...
    textLines,rxtDistType,setDefaultConc,newDefaultRxtConc,...
    setHetConc,hetRxtConc,xllizeParamDefaultCARConc,xllizeHetConc,...
    hetForParams)

% writeHeterogeneityParams.m
% This function writes heterogeneity parameters to the Crystallize3D params
% file.

fid = fopen(heterogeneityParamsOutputFileLocation,'wt');

% Write the first portion of the parameters file
for n = 1:length(textLines(:,1))
    textLine = char(textLines(n,:));
    
    % Write the default reactant amount
    if strfind(textLine,'Default CAR amount')
        if 0 < (newDefaultRxtConc || xllizeParamDefaultCARConc) < 0.000001
            disp('Warning: if reactant concentration is less than 0.000001,')
            disp('Matlab will treat the value as 0.0')
        end
        if setDefaultConc
            line = textLine(1:strfind(textLine,':'));
            fprintf(fid,'%s\t\t\t\t%f\n',line,newDefaultRxtConc);
        else
            line = textLine(1:strfind(textLine,':'));
            fprintf(fid,'%s\t\t\t\t%f\n',line,xllizeParamDefaultCARConc);
        end
        
        % Write the new structures to the parameters file
    elseif strfind(textLine,'structures')
        if 0 < (hetRxtConc || newDefaultRxtConc) < 0.000001
            disp('Warning: if reactant concentration is less than 0.000001,')
            disp('Matlab will treat the value as 0.0')
        end
        
        % Write the line for the new number of structures
        textLineDesc = textLine(1:strfind(textLine,':')+1); % Get all text up to and including ':'
        fprintf(fid,'%s\t\t\t%i\n',char(textLineDesc),char(length(hetForParams(:,1))));
        switch rxtDistType
            case {'uniformBlocks','randomBlocks'}
                %Write the structure lines
                if setHetConc
                    for i = 1:length(hetForParams(:,1))
                        fprintf(fid,'Block (see format codes)					2\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n',...
                            hetForParams(i,1),hetForParams(i,2),hetForParams(i,3),...
                            hetForParams(i,4),hetForParams(i,5),hetForParams(i,6),...
                            hetRxtConc);
                    end
                else
                    for i = 1:length(hetForParams(:,1))
                        fprintf(fid,'Block (see format codes)					2\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n',...
                            hetForParams(i,1),hetForParams(i,2),hetForParams(i,3),...
                            hetForParams(i,4),hetForParams(i,5),hetForParams(i,6),...
                            xllizeHetConc);
                    end
                end
            case {'uniformLayers','randomLayers'}
                %Write the structure lines
                if setHetConc
                    for i = 1:length(hetForParams(:,1))
                        fprintf(fid,'Layer (see format codes)					1\t%i\t%i\t%f\n',...
                            hetForParams(i,1),hetForParams(i,2),hetRxtConc);
                    end
                else
                    for i = 1:length(hetForParams(:,1))
                        fprintf(fid,'Layer (see format codes)					1\t%i\t%i\t%f\n',...
                            hetForParams(i,1),hetForParams(i,2),xllizeHetConc);
                    end
                end
        end
    else
        
        % Write the second portion of the parameters file
        fprintf(fid,'%s\n',char(textLines(n,:)));
    end
end
fclose(fid);

end