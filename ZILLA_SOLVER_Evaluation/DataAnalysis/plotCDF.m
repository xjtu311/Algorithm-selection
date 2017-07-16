function plotCDF(atitle, xtitle, data, mode)
%figure;
%hold on;
xmin=1;
xmax=5000;
ymin=0;
ymax=100;
stepData=100/size(data,1);

if mode==1
semilogx(sort(data), [stepData:stepData:100], 'r-.', 'LineWidth', 3);
end

if mode==2
semilogx(sort(data), [stepData:stepData:100], 'k--', 'LineWidth', 3);
end


if mode==3
semilogx(sort(data), [stepData:stepData:100], 'b-', 'LineWidth', 3);
end


if mode==4
semilogx(sort(data), [stepData:stepData:100], 'g:', 'LineWidth', 3);
end

xlabel (xtitle, 'FontSize', 14);
xlim([xmin xmax]);
ylim([ymin ymax]);
ylabel ('Solved Percentage', 'FontSize', 14);
set(gca,'fontsize',14);
set(gca,'box','on');
% title (atitle, 'FontSize', 18);