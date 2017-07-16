function plotErrorBar(xtitle, ytitle, x, e,color, xlimit, ylimit)
%figure;
%hold on;
if length(xlimit) ==2
    xmin=xlimit(1);
    xmax=xlimit(1);
end

if length(ylimit)==2
    ymin=ylimit(1);
    ymax=ylimit(2);
end
% stepData=100/size(data,1);
%
if color==1
    errorbar(x,e, 'r*');
end

if color==2
    errorbar(x,e, 'bo');
end

if color==3
    errorbar(x,e, 'gs');
end
xlabel (xtitle, 'FontSize', 14);
if length(xlimit) ==2
    xlim([xmin xmax]);
end
if length(ylimit) ==2
    ylim([ymin ymax]);
end
ylabel (ytitle, 'FontSize', 14);
set(gca,'fontsize',14);
set(gca,'box','on');
% title (atitle, 'FontSize', 18);