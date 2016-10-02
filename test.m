
modelArray = zeros(10,8,6);
[modelArray(3:5,2:3,7:9),modelArray(3:5,6:7,4:6),modelArray(5:7,4:5,1:3)] = deal(1)

paddedModelArray = padarray(modelArray,[1 1 1]);

hetIndices = zeros(1,6);
for j = 1:length(paddedModelArray(:,1,1))
    for m = 1:length(paddedModelArray(1,1,:))
        indicesArray = diff(paddedModelArray(j,:,m) ~= 0);
        indicesArrayBeginX = find(indicesArray == 1);
        indicesArrayEndX = find(indicesArray == -1);
        for i = 1:length(indicesArrayBeginX)
            hetIndices = [hetIndices;...
                indicesArrayBeginX(i) j-1 m-1 ...
                indicesArrayEndX(i)-1 j-1 m-1]; % '-1' removes padding
        end
    end
end
hetIndices = hetIndices(2:end,:,:) % get rid of first row of zeros
        
