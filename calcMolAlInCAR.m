function [molAlInCAR...
    ] = calcMolAlInCAR(weightOfCAR,densityOfCAR,numVoxelsOfCAR,voxelVol)
    


molAlInCAR = densityOfCAR * numVoxelsOfCAR * voxelVol / weightOfCAR;
end