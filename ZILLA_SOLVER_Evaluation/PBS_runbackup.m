%% Run Portfolio based algorithm selection for HAL
% input: PBS filename (has all the information)
%        features
% output: solver has been picked
% example: PBS_runbackup([1,2,3], 'ss')
%         Here the set of features are simple features
%% It is a little bit complecatied since we have many different predictor.
function [pick, tElapsed] =PBS_runbackup(f, PBSfilename)
tic; % time information

load(PBSfilename);
tStart = tic;
pick=0;
prediction=0;
if size(f,2)== PBS.NumSimpleFeature
        if length(PBS.FeatModel)>0
            prediction =  applyClassification(PBS.FeatModel, f, 1);
        end
end
if prediction >0
   pick = PBS.Backupsolver;
end
tElapsed = toc(tStart);
end