function transformation = normalData(xTrainAll)

warning off;


constants = determine_constants(xTrainAll, 1);

kept_columns = setdiff(1:size(xTrainAll,2), constants);

xTrainAll = xTrainAll(:,kept_columns);
transformation.kept_columns = kept_columns;

xTrainAllDefault=xTrainAll*0+1;
xTrainAllDefault(find(xTrainAll==-512 | xTrainAll==-1024))=2;


[foo, scale, bias] = determine_transformation(xTrainAll,xTrainAllDefault, 1);
transformation.scale=scale;
transformation.bias=bias;
