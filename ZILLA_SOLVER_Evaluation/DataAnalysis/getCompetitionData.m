dataset='hand';
% dataset='indu';
dataset='rand';
cutoff=5000;
inputfile=sprintf('%s_solver.csv', dataset);
zillafile=sprintf('%s-SATzilla-runtime.csv', dataset);
outputfile=sprintf('%s_all_oracle_zilla_oracle_best.csv', dataset);
solved=3;
runtime=2;
PBS_option.Cutoff=cutoff;

if strcmp(dataset, 'hand')
    allsolver=25;
    take=[2,4,6,7,9,10,12,13,15,18,19,22,23,24,25];
    best=[8];
end
if strcmp(dataset, 'indu')
    allsolver=31;
    take=[1:14,17,18,19, 29];
    best=[11];
end
if strcmp(dataset, 'rand')
    allsolver=17;
    take=[1,2,3,5,6,13,14,15,16];
    best=[9];
end
alldata=csvread(inputfile, 1);
numOutputs=5;
allnames = textread(inputfile, '%s', 1, 'whitespace', '\n', 'bufsize', 10000);
allnames = strread(allnames{1},'%s','whitespace',',');
allnames = deblank(allnames);
namesY = allnames(2:numOutputs:allsolver*5);
allnamesX = allnames(allsolver*5+2:end);
takenamesY=namesY(take);

solvedTable=zeros(size(alldata,1), allsolver);
runtimeTable=alldata(:,runtime:5:5*allsolver);
for i=1:allsolver
    % what we call solved
    foo=solvedTable(:,i)+1;
    foo(alldata(:,i*5-5+solved)==0)=0;
    %     mean(foo)
    foo(alldata(:,i*5-5+runtime)>PBS_option.Cutoff-1)=0;
    alldata(alldata(:,i*5-5+solved)==0,i*5-5+runtime)=PBS_option.Cutoff+1;
    %     mean(foo)
    solvedTable(:,i)=foo;
end
runtimeTable(find(solvedTable==0))=cutoff;
runtimeTable=min(runtimeTable, cutoff);

takeruntimeTable=runtimeTable(:, take);
takesolvedTable=solvedTable(:,take);
solvedbyone=sum(takesolvedTable, 2);
good=find(solvedbyone>0);
length(good)/length(solvedbyone)
takeruntimeTable=takeruntimeTable(good,:);
takesolvedTable=takesolvedTable(good,:);
takesolvedoracle=length(find(sum(takesolvedTable,2)>0));
corrfoo=corr(takeruntimeTable, 'type', 'Spearman');
takeoracle=min(takeruntimeTable, [],2);
oraclecoo=[];
oraclecoo1=[];
takecorr=[];
maxcorr=[];
for i=1:length(corrfoo)
    foo1=corrfoo(:,i);
    foo2=[foo1(1:i-1);foo1(i+1:end)];
    takecorr=[takecorr, foo2];
    maxcorr=[maxcorr, max(foo2)];
    foocoo=min(takeruntimeTable(:,[1:i-1, i+1:end]), [],2);
    foocoo1=length(find(sum(takesolvedTable(:,[1:i-1, i+1:end]),2)>0));
    oraclecoo=[oraclecoo, (mean(foocoo)-mean(takeoracle))/mean(takeoracle)];
    oraclecoo1=[oraclecoo1, (takesolvedoracle-foocoo1)/length(foocoo)];
end

% boxplot(takecorr);
aa=(1-corr(takeruntimeTable, 'type', 'Spearman'))*64;
% IDX=kmeans(takeruntimeTable', 6, 'distance', 'correlation', 'replicates', 1000)
if strcmp(dataset, 'rand')
    order=[3,1,9,2,5,4,6,7,8];
    label={'EagleUP', 'TNM', 'Sparrow',  'Adaptg2wsat11', ...
        'Sattime11', 'MPhaseSAT_M',  'March_rw', 'March_hi', ...
        'Gnovelty+2'};
    zillacoo=[2.2, 0.4, 4.9, -0.4, 0.6, 0.2, 0, 0.2, 0.4 ];
    figure;
    barcoo=[zillacoo', (oraclecoo1*100)'];
    barcoo=barcoo(order,:);
    h=bar(barcoo);
    set(h(2),'facecolor','w')
    xticklabel_rotate([1:9],30,label(order),'interpreter','none');
    ylabel('Cost of Omission (%)', 'fontsize', 16);
    
    %% for random only, I need compute correlation of SAT and UNSAT
    % first find SAT and UNSAT
    satunsatTable=alldata(:,[3:5:5*17]);
    allsolved1=satunsatTable(good, take);
    allsolved=allsolved1;
    allsolved(find(allsolved1==0))=999;
    allsolved=min(allsolved, [],2);
    satid=find(allsolved==1);
    unsatid=find(allsolved==2);
    satruntimeTable=takeruntimeTable(satid,:);
    unsatruntimeTable=takeruntimeTable(unsatid,:);
    satcorr=corr(satruntimeTable, 'type', 'Spearman');
    unsatcorr=corr(unsatruntimeTable, 'type', 'Spearman');
    satcorr=satcorr(order, order);
    unsatcorr=unsatcorr(order, order);
    
    %     figure;
    %     image(aa(order,order), 'CDataMapping', 'direct')
    %     xticklabel_rotate([1:9],30,label(order),'interpreter','none')
    %     colormap(gray)
    %     set(gca,'YTickLabel', label(order));
    %         figure
    %     bar(mean(takesolvedTable(:,order)));
    %      xticklabel_rotate([1:9],30,label(order),'interpreter','none');
    %
    % %      set(gca,'YTickLabel', {''})
    %      ylabel('Solved Percentage', 'fontsize', 16);
    
    meanruntime=mean(takeruntimeTable);
    meansolved=mean(takesolvedTable);
    meancorr=mean(takecorr);
    for i=1:9
        fprintf('%s, %f, %f, %f \n', label{order(i)}, meanruntime(order(i)), meansolved(order(i)), meancorr(order(i)));
    end
    
    
    
end
if strcmp(dataset, 'hand')
    order=  [1     2     3     4     6     5    9  8     10    13     7    11    12    14    15]
    
    label={'MPhaseSAT','Sol', 'QuteRSat', 'CryptoMinisat', 'PicoSAT', ...
        'Glucose2', 'Clasp2', 'Minisat07', 'JMiniSat',  'RestartSAT', 'Clasp1', 'Sathys', 'SApperloT', ...
        'Sattime+', 'Sattime'};
    zillacoo=[3.6, 8.1, 1.4, 0.5, 0.5, 3.6, 2.7, 0, 0.9, 1.4, 1.4 0 1.8 1.4 3.6];
    figure;
    barcoo=[zillacoo', (oraclecoo1*100)'];
    barcoo=barcoo(order,:);
    h=bar(barcoo);
    set(h(2),'facecolor','w')
    xticklabel_rotate([1:15],30,label(order),'interpreter','none');
    ylabel('Cost of Omission (%)', 'fontsize', 16);
    %     figure;
    %     image(aa(order,order), 'CDataMapping', 'direct')
    %     xticklabel_rotate([1:15],30,label(order),'interpreter','none');
    %     colormap(gray)
    %     set(gca, 'YTick', [1:15]);
    %     set(gca,'YTickLabel', label(order));
    %     figure
    %     bar(mean(takesolvedTable(:,order)));
    %      xticklabel_rotate([1:15],30,label(order),'interpreter','none');
    %
    % %      set(gca,'YTickLabel', {''})
    %      ylabel('Solved Percentage', 'fontsize', 16);
    meanruntime=mean(takeruntimeTable);
    meansolved=mean(takesolvedTable);
   meancorr=mean(takecorr);
    for i=1:15
        fprintf('%s, %f, %f, %f \n', label{order(i)}, meanruntime(order(i)), meansolved(order(i)), meancorr(order(i)));
    end
    
end
if strcmp(dataset, 'indu')
    order=[1,4,5,6,9,2,8,3,13,17,7,14,16,10,11,12,15,18];
    
    label={'RestartSAT','Glueminisat', 'Precosat', 'Cirminisat', 'Minisat', ...
        'EBMinisat', 'Contrasat', 'LR GL SHR', 'Minisatagile',  'Glucose1', 'Glucose2', 'EBGlucose', 'Lingeling', ...
        'Minisat_psm', 'CryptoMinisat', 'Rcl', 'MPhaseSAT64', 'QuteRSat'};
    zillacoo=[0.4, 0.8, 0.4, 0.0, -0.4, 0, 1.2, 0.8, 0.4, 0, 0.4, 0.8, 0.8, 1.2 0.4 0.4 1.6 0.8];
    figure;
    barcoo=[zillacoo', (oraclecoo1*100)'];
    barcoo=barcoo(order,:);
    h=bar(barcoo);
    set(h(2),'facecolor','w')
    xticklabel_rotate([1:18],30,label(order),'interpreter','none');
    ylabel('Cost of Omission (%)', 'fontsize', 16);
    %     figure;
    %     image(aa(order,order), 'CDataMapping', 'direct')
    %     xticklabel_rotate([1:18],30, label(order),'interpreter','none');
    %     colormap(gray)
    %     set(gca, 'YTick', [1:18]);
    %     set(gca,'YTickLabel', label(order));
    %         figure
    %     bar(mean(takesolvedTable(:,order)));
    %      xticklabel_rotate([1:18],30,label(order),'interpreter','none');
    % %      set(gca,'YTickLabel', {''})
    %      ylabel('Solved Percentage', 'fontsize', 16);
    meanruntime=mean(takeruntimeTable);
    meansolved=mean(takesolvedTable);
    meancorr=mean(takecorr);
    for i=1:18
        fprintf('%s, %f, %f, %f \n', label{order(i)}, meanruntime(order(i)), meansolved(order(i)), meancorr(order(i)));
    end
    
end
% set(gca,'XTickLabel', {''})
% ylim([-0.2,1.05]);
% % set(findobj(gca,'Type','text'),'FontSize',12)
% ylabel('Correlation Coefficient', 'fontsize', 16);

% figure
% bar(oraclecoo1);
% ylabel('Cost of Omission (Oracle)', 'fontsize', 18);
% figure
% plot(maxcorr, oraclecoo, 'r.');



%% plot satzilla's performance compare to other solvers
% o1=min(runtimeTable,[],2);
%
% o2=min(runtimeTable(:,take),[],2);
% alloutput=[o1,o2, runtimeTable(:,best)];
%
% zillatime=csvread(zillafile);
% alloutput=min(alloutput, cutoff+1);
% zillatime=min(zillatime, cutoff+1);
%
% csvwrite(outputfile,[alloutput, zillatime]);
%
% alloutput=max(alloutput, 0.001);
% zillatime=max(zillatime,0.001);
%
% if strcmp(dataset, 'hand') || strcmp(dataset, 'rand')
% figure;
% plotRuntime('SATzilla(HAND) Runtime (CPU seconds)', 'Virtual Best Solver Runtime (CPU seconds)', [zillatime, alloutput(:,1)]);
% figure;
% plotRuntime('SATzilla(HAND) Runtime (CPU seconds)', 'Oralce(HAND) Runtime (CPU seconds)', [zillatime, alloutput(:,2)]);
% figure;
% plotRuntime('SATzilla(HAND) Runtime (CPU seconds)', '3S Runtime (CPU seconds)', [zillatime, alloutput(:,3)]);
% figure;
% alloutput(alloutput>cutoff)=cutoff*100;
% zillatime(zillatime>cutoff)=cutoff*100;
%
% plotCDF('', 'Runtime (CPU seconds)', alloutput(:,1), 1);
% hold on;
% plotCDF('', 'Runtime (CPU seconds)', alloutput(:,2), 2);
% plotCDF('', 'Runtime (CPU seconds)', zillatime, 3);
% plotCDF('', 'Runtime (CPU seconds)', alloutput(:,3), 4);
% hold off;
% end
% if strcmp(dataset, 'indu')
% figure;
% plotRuntime('SATzilla(INDU) Runtime (CPU seconds)', 'Virtual Best Solver Runtime (CPU seconds)', [zillatime, alloutput(:,1)]);
% figure;
% plotRuntime('SATzilla(INDU) Runtime (CPU seconds)', 'Oralce(INDU) Runtime (CPU seconds)', [zillatime, alloutput(:,2)]);
% figure;
% plotRuntime('SATzilla(INDU) Runtime (CPU seconds)', 'GLUCOSE 2 Runtime (CPU seconds)', [zillatime, alloutput(:,3)]);
% figure;
% alloutput(alloutput>cutoff)=cutoff*100;
% zillatime(zillatime>cutoff)=cutoff*100;
%
% plotCDF('', 'Runtime (CPU seconds)', alloutput(:,1), 1);
% hold on;
% plotCDF('', 'Runtime (CPU seconds)', alloutput(:,2), 2);
% plotCDF('', 'Runtime (CPU seconds)', zillatime, 3);
% plotCDF('', 'Runtime (CPU seconds)', alloutput(:,3), 4);
% hold off;
% end