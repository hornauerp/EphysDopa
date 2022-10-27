function accuracy_barplot(timeline_accuracy)
capsz = 3;
fontsz = 8;
bar([mean(timeline_accuracy(:,1:5)) 0 mean(timeline_accuracy(:,6))],'k','FaceAlpha',0.8)
hold on
err = [std(timeline_accuracy(:,1:5)) 0 std(timeline_accuracy(:,6))];
er = errorbar(1:7,[mean(timeline_accuracy(:,1:5)) 0 mean(timeline_accuracy(:,6))],err,err,'k','linestyle','none','CapSize',capsz);
% yl = yline(chance_line,'k:','LineWidth',2,'Alpha',0.3);
% annotation('line',[0.245 0.25],[0.04 0.06])
% annotation('line',[0.24 0.245],[0.04 0.06])
xlabel('Week')
set(gca,'XTickLabel',{1:5,'All'})
set(gca,'XTick',[1:5,7])
set(gca,'YLim',[0.5 1.05])
set(gca,'YTick',[0.4 0.6 0.8 1])
ylabel('Accuracy')
box off
set(gca,'FontSize',fontsz)