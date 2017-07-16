%% Run Portfolio based algorithm selection for HAL
% input: PBS filename (has all the information)
%        features
% output: solver has been picked
% example: PBS_run([1,2,3], 'ss')
%
%% It is a little bit complecatied since we have many different predictor.
function [pick, tElapsed] =PBS_run(f, PBSfilename)
tic; % time information
tStart = tic;

load(PBSfilename);
pick = PBS.Backupsolver;
if size(f,2)== PBS.NumFeature
    if length(PBS.SolverPicked)>1
        if strcmp(PBS.Predictor, 'LR') || strcmp(PBS.Predictor, 'RF')
            % make prediction for all solvers
            prediction=[];
            for i =1:length(PBS.SolverPicked)
                prediction(i) = applyModel(PBS.Model{PBS.SolverPicked(i)}, f, 0, 0, 0);
            end
            [a,b]=sort(prediction,'descend');
            pick=PBS.SolverPicked(b);
        end
        if strcmp(PBS.Predictor,'SVM') || strcmp(PBS.Predictor, 'DT') || strcmp(PBS.Predictor, 'DF')
            % make prediction for all pairwise solvers
            prediction=[];
            predid=1;
            for ii = 1:length(PBS.SolverPicked)
                for jj =ii+1:length(PBS.SolverPicked)
                    modelid=find(PBS.PredictItem1==PBS.SolverPicked(ii) &  PBS.PredictItem2==PBS.SolverPicked(jj));
                    prediction(predid) = applyClassification(PBS.Model{modelid}, f, 0);
                    predid=predid+1;
                end
            end
            % then, we need find out which one is the best
            tmppick = selectSolverSVM(prediction, length(PBS.SolverPicked), 1);
            pick = PBS.SolverPicked(tmppick);
        end
    else
        pick=PBS.SolverPicked;
    end
end
tElapsed = toc(tStart);

end