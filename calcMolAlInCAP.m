function [molAlInCAP...
    ] = calcMolAlInCAP(AlInCAP,weightOfCAP,densityOfCAP,modeOfCAP)
    


molAlInCAP = AlInCAP * densityOfCAP * modeOfCAP / weightOfCAP;
end