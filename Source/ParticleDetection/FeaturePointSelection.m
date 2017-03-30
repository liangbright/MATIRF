function [SelectedFeaturePoint SelectedIdxList FeatureClassifier]=FeaturePointSelection(FeaturePoint, Feature_SNR_Thresh)

if isempty(FeaturePoint)
    SelectedFeaturePoint=[];
    SelectedIdxList=[];
    return;
end


Feature=FeaturePoint([3, 4], :);

FeatureClassifier=GMMClassifierClass();
Label=FeatureClassifier.Build(Feature, Feature_SNR_Thresh);

SelectedIdxList=find(Label>0);
SelectedFeaturePoint=FeaturePoint(:, SelectedIdxList);

end