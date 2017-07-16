function plotCDF(xtitle, ytitle, data)
%figure;
%hold on;
xmin=0.001;
xmax=5100;
ymin=0.001;
ymax=5100;

loglog(data(:,1), data(:,2), 'r.', 'MarkerSize', 10);

xlabel (xtitle, 'FontSize', 14);
xlim([xmin xmax]);
ylim([ymin ymax]);
ylabel (ytitle, 'FontSize', 14);
set(gca,'fontsize',14);
set(gca,'box','on');
hold on;
line([xmin,xmax], [ymin, ymax],'LineWidth', 4 );
hold off;

% title (atitle, 'FontSize', 18);