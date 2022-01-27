function heatmap_plot(sc_timeline_data,sc_features,nw_timeline_data,nw_features,idx1,idx2,leg)
heatmap_data = cat(3,sc_timeline_data,nw_timeline_data);

heatmap_data(:,:,strcmp(sc_features,'PACF')) = heatmap_data(:,:,strcmp(sc_features,'PACF'))+1; %Make PACF values only positive
heatmap_log_ratio = log(squeeze(heatmap_data(idx2,:,:)./heatmap_data(idx1,:,:)));
heatmap_log_ratio(heatmap_log_ratio>log(3)) = log(3); %Limit displayed value range
heatmap_log_ratio(heatmap_log_ratio<log(1/3)) = log(1/3);
heatmap_log_ratio = [heatmap_log_ratio(:,1:length(sc_features)), zeros(5,1),heatmap_log_ratio(:,length(sc_features)+1:end)];
labels = string([sc_features;"";nw_features]);

boxsize = 0.02; 
diff_colorbar = 0.01;
cb_height = 0.01;
fontsz = 7;
cm = othercolor('RdBu9',100); cm((size(cm,1)/2)-2:size(cm,1)/2+2,:) = [];

figure('Color','w','Position',[1200 100 300,600]);
%Heatmap subplot
s_hmap = axes;
hmap = imagesc(heatmap_log_ratio');
hold on
%Mark plotted features in heatmap
s_hmap.TickLength = [0 0];
s_hmap.FontSize = fontsz;
s_hmap.YTick = 1:length(labels);
s_hmap.YTickLabel = labels;
s_hmap.XTick = 1:5;
s_hmap.XLabel.String = "Week";
colormap(cm)
cb = colorbar('Position',[s_hmap.Position(1) s_hmap.Position(2)+s_hmap.Position(4)+diff_colorbar, s_hmap.Position(3), cb_height],...
    'Location','northoutside','Limits',[log(1/3),log(3)]);
cb.Ticks = [log(1/3),0,log(3)];
cb.TickLabels = ["<0.3", "1", ">3"];
cb.Ruler.TickLabelGapOffset = -1;
if any(size(leg)==1)
    ttl_str = upper(strjoin({leg{idx2},'/',leg{idx1}},''));
else
    ttl_str = upper(strjoin({leg{2,idx2},' ',leg{1,idx2},'/',leg{2,idx1},' ',leg{1,idx1}},''));
end
cb.Title.String = ['\bf' ttl_str];
% cb.Title.Position = [21 14 0];
feature_sep1 = yline(s_hmap,length(sc_features)+0.5,'LineWidth',1);
feature_sep2 = yline(s_hmap,length(sc_features)+1.5,'LineWidth',1);