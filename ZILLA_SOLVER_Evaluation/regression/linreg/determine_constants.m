function constants = determine_constants(xTrain, trafo)

xTrainDefault=xTrain*0+1;
xTrainDefault(find(xTrain==-512 | xTrain==-1024))=2;

nonconstant = determine_transformation(xTrain, xTrainDefault, trafo);
constants = setdiff(1:size(xTrain,2), nonconstant);

fprintf('Discarding %i constant features of %i in total.\n', length(constants), size(xTrain,2));