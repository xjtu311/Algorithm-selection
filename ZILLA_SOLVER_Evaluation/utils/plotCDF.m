function plotCDF(atitle, xtitle, data)
%figure;
%hold on;
xmin=1;
xmax=4000;
ymin=0;
ymax=100;
stepData=100/size(data,1);

semilogx(sort(data), [stepData:stepData:100], 'r-');
xlabel (xtitle, 'FontSize', 14);
xlim([xmin xmax]);
ylim([ymin ymax]);
ylabel ('Solved Percentage', 'FontSize', 14);
set(gca,'fontsize',14);
set(gca,'box','on');
% title (atitle, 'FontSize', 18);