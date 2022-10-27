function f = plot_features(sel_vars,data_matrix,data_sd,tps,leg)

n_groups = size(data_matrix,1);
c = othercolor('RdBu4',n_groups);

if n_groups == 2
    c([1 2],:) = c;
elseif n_groups == 4
    c([1 2 4 3],:) = c;
elseif n_groups == 6
    c([1 2 3 6 5 4],:)= c([1:6],:);
end

if any(size(leg)==1)
    leg_entry = upper(leg);
else
    
    for i = 1:size(leg,2)
        leg_entry(i) = string(strjoin({upper(leg{2,i}),leg{1,i}}));
    end
end

f = figure('Name','Single-cell Features','Color','w','Position',[100 100 1500 1000]);
tiledlayout('flow','TileSpacing','tight')
% c([1 3 2 6 4 5],:) = c;

jitter = linspace(-1,1,n_groups);
fontsz = 8;
for v = 1:length(sel_vars)
    tl_plot{v} = nexttile;
    
    for i = 1:n_groups
        x = tps+1+jitter(i);
        y = data_matrix(i,1:length(x),v); %Quick fix, has to be corrected to show the selected time points
        xx = linspace(min(x),max(x),100); yy = pchip(x,y,xx);
        plot(xx,yy,'Color',c(i,:),'HandleVisibility','off')
        hold on
        errorbar(tps+1+jitter(i),data_matrix(i,1:length(x),v),data_sd(i,1:length(x),v),...
            'LineWidth',1,'Color',c(i,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
            'MarkerFaceColor',c(i,:));
        set(gca,'FontSize',fontsz)
        marg = get(gca,'ylabel');
        set(marg,'Margin',3)
    end
    title(sel_vars{v})
    xlabel('DIV')
%     ylabel('Hz')
    box off
end
l = legend(leg_entry,'Box','off');
l.Position = [0.02 0.8 0.1 0.1];
sgtitle(sprintf('Development until Day %d',tps(end)+1))