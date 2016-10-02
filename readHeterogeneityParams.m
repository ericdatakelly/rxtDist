function [rxtDistType,heterogeneityLengthX,heterogeneityLengthY,...
    heterogeneityLengthZ,hetRandomSeed,fixedGapX,fixedGapY,fixedGapZ,...
    setDefaultConc,newDefaultRxtConc,setHetConc,hetRxtConc...
    ] = readHeterogeneityParams(name)

% readHeterogeneityParams.m
% This function reads heterogeneity parameters for input to Crystallize3D


fid = fopen(name);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':')+1:length(line)); % get everything after ':' to end of line
rxtDistType = sscanf(line,'%s',1); % get the string of characters
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
heterogeneityLengthX = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
heterogeneityLengthY = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
heterogeneityLengthZ = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
hetRandomSeed = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
fixedGapX = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
fixedGapY = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
fixedGapZ = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
setDefaultConc = sscanf(line,'%*c %i',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
newDefaultRxtConc = sscanf(line,'%*c %g',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
setHetConc = sscanf(line,'%*c %i',1);
line = fgetl(fid);
line = line(strfind(line,':'):length(line));
hetRxtConc = sscanf(line,'%*c %g',1);
fclose(fid);
