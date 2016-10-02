function [] = rxtDist(xllizeParamsName,inputHetParamsName,...
    outputXllizeParamsName)

% rxtDist.m
% syntax: rxtDist('filename1','filename2','filename3')
% Example: rxtDist('params_PM3_140717_1437.txt','input.txt','output.txt')
% filename1: Crystallize parameter file name
% filename2: Heterogeneity parameter file name
% filename3: Name of new parameters file written with new heterogeneity
% parameters
% File locations: filename1 and filename2 must be within the Params folder 
% inside the current directory.  The output file will be written to the
% same directory.
%
% This program creates various types of reactant distributions for use in
% Crystallize3D (Ketcham and Carlson, 2012)
% More description...


% Set up file names and paths
paramsPath = strcat(cd,'\Params\');
heterogeneityParamsFileLocation = strcat(paramsPath,inputHetParamsName);
heterogeneityParamsOutputFileLocation = strcat(paramsPath,outputXllizeParamsName);
crystallizeParamsInputFileLocation = strcat(paramsPath,xllizeParamsName);

fid = fopen(crystallizeParamsInputFileLocation);
if fid == -1
    disp('---------------------------------------------');
    disp('   Error in rxtDist.m');
    disp('   Cannot find Crystallize3D parameters file');
    disp(['   "',xllizeParamsName,'" is not in the directory']);
    disp('---------------------------------------------');
    return
else
    fclose(fid);
end


% Read Crystallize3D parameters file
[textLines,rxtStructures,xllizeParamStructureType,parameterArray...
    ] = readXllizeParams(crystallizeParamsInputFileLocation);
    
% Read the heterogeneity parameters
[rxtDistType,heterogeneityLengthX,heterogeneityLengthY,...
    heterogeneityLengthZ,hetRandomSeed,gapX,gapY,gapZ,...
    setDefaultConc,newDefaultRxtConc,setHetConc,hetRxtConc...
    ] = readHeterogeneityParams(heterogeneityParamsFileLocation);


% Establish several parmeters
% currently the headers take up lines in the array as zeros
voxelEdgeLength = parameterArray(20);
xllizeParamMaxX = parameterArray(22);
xllizeParamMaxY = parameterArray(23);
xllizeParamMaxZ = parameterArray(24);
xllizeParamDimensions = parameterArray(25);
xllizeParamDefaultRxtConc = parameterArray(26);
xllizeParamNumStructures = parameterArray(28);
voxelVol = voxelEdgeLength^3;


% Determine amount of reactant in Crystallize3D model
[xllizeRxtAmount,xllizeHetConc...
    ] = calcXllizeRxtAmt(voxelVol,xllizeParamMaxX,xllizeParamMaxY,...
    xllizeParamMaxZ,xllizeParamDimensions,xllizeParamDefaultRxtConc,...
    xllizeParamStructureType,xllizeParamNumStructures,rxtStructures);


% Check that some parameters are reasonable
if (heterogeneityLengthX > xllizeParamMaxX - 1 || ...
        heterogeneityLengthY > xllizeParamMaxY - 1 || ...
        heterogeneityLengthZ > xllizeParamMaxZ - 1)
    disp('Heterogeneity dimension must be less than the model dimension.')
    disp('Adjust dimension and try again.')
    return
end
if heterogeneityLengthX == 0 || ...
        heterogeneityLengthY == 0 || ...
        heterogeneityLengthZ == 0
    disp('Heterogeneity dimension must be an integer greater than zero.')
    disp('Adjust dimension and try again.')
    return
end
if gapX > xllizeParamMaxX || ...
        gapY > xllizeParamMaxY || ...
        gapZ > xllizeParamMaxZ
    disp('Gap length must be less than the model dimension.')
    disp('Adjust dimension and try again.')
    return
end


% Calculate the new heterogeneities
switch rxtDistType
    case 'uniformLayers'
        [hetForParams...
            ] = calcLayersUniform(xllizeParamMaxX,xllizeParamMaxY,...
            xllizeParamMaxZ,gapX,heterogeneityLengthX,...
            xllizeParamDimensions);
    case 'randomLayers'
        [hetForParams...
            ] = calcLayersRandom(xllizeParamMaxX,xllizeParamMaxY,...
            xllizeParamMaxZ,heterogeneityLengthX,xllizeHetConc,...
            xllizeParamDefaultRxtConc,voxelVol,hetRandomSeed,...
            xllizeRxtAmount,xllizeParamDimensions);
    case 'uniformBlocks'
        [hetForParams...
            ] = calcBlocksUniform(xllizeParamMaxX,xllizeParamMaxY,...
            xllizeParamMaxZ,gapX,gapY,gapZ,heterogeneityLengthX,...
            heterogeneityLengthY,heterogeneityLengthZ,xllizeParamDimensions);
    case 'randomBlocks'
        [hetForParams...
            ] = calcBlocksRandom(xllizeParamMaxX,xllizeParamMaxY,...
            xllizeParamMaxZ,heterogeneityLengthX,heterogeneityLengthY,...
            heterogeneityLengthZ,xllizeHetConc,xllizeParamDefaultRxtConc,...
            voxelVol,hetRandomSeed,xllizeRxtAmount,xllizeParamDimensions);
end
hetForParams

% Write the coordinates of the heterogenieities to the params file with the
% previously read xllize params (porosity, filename, etc.).
writeHeterogeneityParams(heterogeneityParamsOutputFileLocation,...
    textLines,rxtDistType,setDefaultConc,newDefaultRxtConc,...
    setHetConc,hetRxtConc,xllizeParamDefaultRxtConc,xllizeHetConc,...
    hetForParams);

