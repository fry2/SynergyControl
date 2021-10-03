varHolder = whos;

if sum(contains({varHolder.name},{'BackMean', 'BackRaw'})) == 0
    [Stride, StrideContact, NewTime, BackMean,BackRaw] = CompareAnimalAndSim(0);
end

extrapolatedData = zeros(size(BackRaw));
baseTable = array2table(BackRaw(:,:,i));

for i = 1:3    
    baseTable.Var9 = BackMean(:,i);
    [trainedModel, validationRMSE] = trainRegressionModel(baseTable);
    for j = 1:8
        if any(isnan(BackRaw(:,j,i)))
            baseTable.Var9 = BackRaw(:,j,i);
            extrapolatedData(:,j,i) = trainedModel.predictFcn(baseTable);
        else
            extrapolatedData(:,j,i) = BackRaw(:,j,i);
        end
    end
end