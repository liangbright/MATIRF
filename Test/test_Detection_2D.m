Iavg=sum(ImageStack, 3)/AngleNum;
TypicalRadius=PSFSigma;
FeaturePoint_Raw=FeaturePointDetection(Iavg, CellLabel, TypicalRadius);
FeaturePoint=FeaturePointSelection(FeaturePoint_Raw, 2);
x_Init=FeaturePoint(1,:);
y_Init=FeaturePoint(2,:);