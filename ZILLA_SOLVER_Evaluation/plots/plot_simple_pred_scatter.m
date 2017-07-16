function plot_simple_pred_scatter(y, y_cross, y_cross_var, cens, rmse, cc, ll, figure_prefix, title_prefix, logModel)
nCross = length(y_cross);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% observed values vs. cross-validated prediction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cens_idx=find(cens==1);
uncens_idx=find(cens==0);

%sfigure(6);
sfigure;
yplot = y;
if logModel
    ycross_plot = 10.^y_cross;
else
    ycross_plot = y_cross;
    mini = min(min(y), min(y_cross))/1.1;
    maxi = max(max(y), max(y_cross))*1.1;
end
mini = min(min(y), min(ycross_plot))/1.1;
maxi = max(max(y), max(ycross_plot))*1.1;
    
hold off
hE     = loglog(yplot, ycross_plot);
hold on
hXLabel = xlabel('true runtime');
hYLabel = ylabel('predicted runtime');
hLine = line([mini, maxi],[mini,maxi]);

set(hLine                         , ...
  'Color'           , [0 0 .5]    , ...
  'LineWidth'       , 2           );

set([hE]                     , ...
  'LineStyle'       , 'none'      , ...
  'Color'           , [.3 .3 .3]  , ...
  'LineWidth'       , 1           , ...
  'MarkerSize'      , 3           , ...
  'Marker'          , 'o'         , ...
  'MarkerFaceColor' , [.7 .7 .7]  , ...
  'MarkerEdgeColor' , [.2 .2 .2]);

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'FontSize'    , 14, ...
  'XTick'       , 10.^[-10:10], ...
  'YTick'       , 10.^[-10:10], ...
  'LineWidth'   , 2         );
%  'TickLength'  , [.02 .02] , ...
%  'YGrid'       , 'on'      , ...
%  'YTick'       , 0:500:2500, ...

%Xcolor [.3 .3 .3]

% if logModel
%     set(gca, ...
% end

axis([mini, maxi, mini, maxi]);

set([hXLabel, hYLabel]  , ...
    'FontSize'   , 18          );
set(gcf, 'Outerposition', [0,0,500,500]);

if ~strcmp(figure_prefix, '')
    set(gcf, 'PaperPositionMode', 'auto');
%    filename = strcat(figure_prefix, 'pred.eps');
%    fprintf(strcat(['Saving plot to ', filename]));
%    print('-depsc2', filename);
    filename = strcat(figure_prefix, '-pred.pdf');
    if exist(filename, 'file')
        delete(filename)
    end

%    export_fig(filename);

    saveas(gcf, strcat(figure_prefix,'-pred.fig'));
%    close;
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % standardized residual plot from the EGO paper.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h=sfigure(7);
% hold off
% stand_res = 1e10 * ones(length(y_cross_var),1);
% nz_idx = find(y_cross_var>0);
% stand_res(nz_idx) = (yplot(nz_idx)-y_cross(nz_idx))./(y_cross_var(nz_idx).^0.5);
% 
% % if length(y) > 100
% %     idx = 1:100;
% %     uncens_idx_sub = uncens_idx(intersect(uncens_idx,idx));
% %     cens_idx_sub = cens_idx(intersect(cens_idx,idx));
% %     hUncens = semilogx(y(uncens_idx_sub), stand_res(uncens_idx_sub), 'ro');
% %     hold on
% %     hCens = semilogx(y(cens_idx_sub), stand_res(cens_idx_sub), 'rx');
% %     
% %     idx = 101:length(y);
% %     uncens_idx_sub = uncens_idx(intersect(uncens_idx,idx));
% %     cens_idx_sub = cens_idx(intersect(cens_idx,idx))-100;
% %     hUncens = semilogx(y(uncens_idx_sub), stand_res(uncens_idx_sub), 'go');
% %     hold on
% %     hCens = semilogx(y(cens_idx_sub), stand_res(cens_idx_sub), 'gx');
% %     
% %     hXLabel = xlabel('observed value');
% %     hYLabel = ylabel('standardized residual cross-validation error')
% % else
%     hUncens = semilogx(y(uncens_idx), stand_res(uncens_idx), 'ko');
%     hold on
%     hCens = semilogx(y(cens_idx), stand_res(cens_idx), 'kx');
%     if meanRT
%         %hXLabel = xlabel('true penalized average runtime [s]');
%         hXLabel = xlabel('true mean response');
%         hYLabel = ylabel('standardized residual error');
%     else
%         %hXLabel = xlabel('observed penalized runtime [s]');
%         hXLabel = xlabel('observed response');
%         hYLabel = ylabel('standardized residual CV error');
%     end
% % end
% hLine1 = line([min(y)/2, max(y)*2],[-3,-3]);
% hLine2 = line([min(y)/2, max(y)*2],[ 3, 3]);
% 
% set(hUncens                       , ...
%   'LineStyle'       , 'none'      , ...
%   'Marker'          , 'o'         , ...
%   'MarkerSize'      , 6           , ...
%   'MarkerFaceColor' , [1 1 1]); %[.7 .7 .7]  );
% %  'MarkerEdgeColor' , [.2 .2 .2]  , ...
% 
% set([hLine1, hLine2]              , ...
%   'Color'           , [0 0 .5]    , ...
%   'LineWidth'       , 2           );
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'XTick'       , 10.^[-10:10], ...
%       'FontSize'   , 14, ...  
%   'LineWidth'   , 1         );
% %  'TickLength'  , [.02 .02] , ...
% %  'YTick'       , 0:500:2500, ...
% 
% set(gcf, 'Outerposition', [0,0,500,500]);
% %  'XColor'      , [.3 .3 .3], ...
% %  'YColor'      , [.3 .3 .3], ...
% 
% 
% axis([min(y)/2, max(y)*2, min(min(stand_res)-0.01, -3.5), max(max(stand_res)+0.01, 3.5)]);
% 
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 18          );
% 
% if ~strcmp(figure_prefix, '')
%     set(gcf, 'PaperPositionMode', 'auto');
% %    filename = strcat(figure_prefix, 'err.eps');
% %    fprintf(strcat(['Saving plot to ', filename]));
% %    print('-depsc2', filename);
%     filename = strcat(figure_prefix, 'err.pdf');
%     if exist(filename, 'file')
%         delete(filename)
%     end
%     export_fig(filename);
%     saveas(gcf, strcat(figure_prefix,'err.fig'))
% %    close;
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % standardized normal quantile plot from the EGO paper.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nrand = randn(1,1000);
% real_q = [];
% my_q = [];
% %  stand_res = randn(1,40);
% for i=1:length(uncens_idx)
%     real_q(i) = quantile(nrand, (i+0.0)/(nCross+1.0));
%     my_q(i) = quantile(stand_res(uncens_idx), (i+0.0)/(nCross+1.0));
% end
% h=sfigure(8);
% hold off
% hDots = plot(real_q, my_q, 'ko', 'MarkerSize', 6);
% 
% hold on
% hXLabel = xlabel('standard normal quantile');
% hYLabel = ylabel('standardized residual quantile');
% mini = min(min(real_q), min(my_q)) - 0.5;
% maxi = max(max(real_q), max(my_q)) + 0.5;
% hLine = line([mini, maxi],[mini,maxi]);
% 
% set(hLine           , ...
%   'Color'           , [0 0 .5]    , ...
%   'LineWidth'       , 2           );
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'XGrid'       , 'on'      , ...
%       'FontSize'   , 14, ...  
%   'LineWidth'   , 1         );
% %  'TickLength'  , [.02 .02] , ...
% %  'YTick'       , 0:500:2500, ...
% 
% %  'XColor'      , [.3 .3 .3], ...
% %  'YColor'      , [.3 .3 .3], ...
% 
% 
% axis([mini, maxi, mini, maxi]);
% set(gcf, 'Outerposition', [0,0,500,500]);
%     
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 18          );
% 
% if ~strcmp(figure_prefix, '')
%     set(gcf, 'PaperPositionMode', 'auto');
% %    filename = strcat(figure_prefix, 'qq.eps');
% %    fprintf(strcat(['Saving plot to ', filename]));
% %    print('-depsc2', filename);
%     filename = strcat(figure_prefix, 'qq.pdf');
%     if exist(filename, 'file')
%         delete(filename)
%     end
%     export_fig(filename);
% 
%     saveas(gcf, strcat(figure_prefix,'qq.fig'));
% %    close;
% end