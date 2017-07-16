function prediction = applyClassification(model, X, isclean)
global svmpath;
if ~isclean
    if ~isfield(model, 'transformation')
        error('ERROR. The model was learnt on clean data, it doesnt have a transformation for unclean data.');
    end
    X = formatnormalData(X, model.transformation);
end

switch model.type
    case 'SVM'
        testfile=sprintf('%stest%f', svmpath, rand());
        resultfile=sprintf('%sresult%f', svmpath, rand());
        fout = fopen(testfile, 'w');
        for uu=1:size(X,1)
            fprintf(fout, '%d cost:%f ', 1,  1.0);
            for ww=1:size(X,2)
                fprintf(fout, '%d:%f ', ww, X(uu, ww));
            end
            fprintf(fout, '\n');

        end
        fclose(fout);
        mycmd=sprintf('!%ssvm_classify %s %s %s', svmpath, testfile, model.modelfile, resultfile);
        eval(mycmd);
        prediction=csvread(resultfile);
        delete(testfile);
        delete(resultfile);
    case 'DT'
        predclass=eval(model.tree, X);
        prediction=str2num(char(predclass));
    case 'DF'
        pred_test=[];
        for num_tree=1:length(model.trees)
            pred_test(:, num_tree) = str2num(char(eval(model.trees{num_tree}, X)));
        end
        prediction = mean(pred_test,2);

        prediction(find(prediction>0))=1;
        prediction(find(prediction<0))=-1;
end
end