function xData = formatnormalData(xData, transformation)

%=== Remove constant features.
xData = xData(:,transformation.kept_columns);
%=== deal with default feature
xDataDefault=xData*0+1;
xDataDefault(find(xData==-512 | xData==-1024))=2;
%=== Normalize data.
%%% I just turn this off for mixture of expert %%%%
xData = (xData - repmat(transformation.bias, [size(xData,1),1])) ./ repmat(transformation.scale, [size(xData,1),1]);
xData(find(xDataDefault>1))=0;

