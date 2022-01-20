function [silhouette_coef,sum_explained] = plotPCA(train_mat,train_color,test_mat,test_color,viewAngle,markerSize,onlyArea)
arguments
train_mat {mustBeNumeric}
train_color (1,:)
test_mat {mustBeNumeric} = [];
test_color (1,:) = [];
viewAngle = [];
markerSize (1,1) {mustBeNumeric} = 20;
onlyArea (1,1) = 0;
end
%% Calculate PC space
mu = mean(train_mat,'omitnan');
sd = std(train_mat,'omitnan');
norm_train = (train_mat-mu)./sd;

[coeff,train_score,~,~,explained] = pca(norm_train);
if ~isempty(test_mat)
norm_treat = (test_mat-mu)./sd;
treatment_score = (norm_treat-mu)*coeff;
else
    treatment_score = [];
end
pcs = [1:3];
score = [train_score; treatment_score];
score(isnan(score)) = 0;
color_idx = [train_color test_color];
n_clust = unique(color_idx);

cb = othercolor('RdBu4',max(color_idx));
if length(n_clust) == 2
    cb([3 2 1],:) = cb;
elseif length(n_clust) == 4
cb([3,4,2,1],:) = cb;
end
centroids = nan(length(n_clust),3);
x_plot = cell(1,4); y_plot = cell(1,4); z_plot = cell(1,4);
x_sel = cell(1,4); y_sel = cell(1,4); z_sel = cell(1,4);
for i = n_clust
    x_plot{i} = score(color_idx==i,pcs(1));
    y_plot{i} = score(color_idx==i,pcs(2));
    z_plot{i} = score(color_idx==i,pcs(3));
        centroids(i,:) = [mean(x_plot{i}),mean(y_plot{i}),mean(z_plot{i})];
        dists = pdist2([x_plot{i} y_plot{i} z_plot{i}],centroids(i,:));
    if length(x_plot{i})>4
        outlier_idx = isoutlier(dists,'ThresholdFactor',10);
    else
        outlier_idx = zeros(1,length(x_plot{i}));
    end
    x_sel{i} = x_plot{i}(~outlier_idx);
    y_sel{i} = y_plot{i}(~outlier_idx);
    z_sel{i} = z_plot{i}(~outlier_idx);
end
for i = n_clust
    if length(n_clust) == 4 && ismember(i,[2 4]) || length(n_clust)==2 || size(score,1)<20
    scatter3(x_plot{i},y_plot{i},z_plot{i},markerSize,cb(i,:),'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.8)
    end
    hold on
%     s = scatter3(mean(x_plot),mean(y_plot),mean(z_plot),50,cb(i,:),'filled','MarkerEdgeColor','k','LineWidth',1);
%     set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
if length(x_sel{i})>3
    x = x_sel{i}; y = y_sel{i};z = z_sel{i};
    j = boundary(x,y,z,0);
%     s = interparc(1000,x([1:end 1]),y([1:end 1]),z([1:end 1]),'csape');
%     j = boundary(s(:,1),s(:,2),s(:,3),0);
%     t = trisurf(j,s(:,1),s(:,2),s(:,3), 'FaceColor',cb(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.0);
    t = trisurf(j,x,y,z, 'FaceColor',cb(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.1);
else
%     t = fill3(x_plot,y_plot,z_plot,'FaceColor',cb(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.2);
    t = patch('Faces',1:3,'Vertices',[x_plot{i} y_plot{i} z_plot{i}], 'FaceColor',cb(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.2);
end
if ~onlyArea || (length(n_clust) == 4 && ~ismember(i,[1 3]))
    set(get(get(t,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
end
% p = plot3(centroids(n_clust,1),centroids(n_clust,2),centroids(n_clust,3),'k--');
% p.Color(4) = 0.5;
xlabel(sprintf('PC %d',pcs(1)))
ylabel(sprintf('PC %d',pcs(2)))
zlabel(sprintf('PC %d',pcs(3)))
axis tight
% if length(n_clust) == 2
% legend({'WT','PD'})
% elseif length(n_clust) == 4
% legend({'WT','WT+LNA','PD','PD+LNA'})
% end
if ~isempty(viewAngle)
view(viewAngle)
end

silhouette_coef = mean(silhouette(score(:,1:3),color_idx,'Euclidean'));
sum_explained = sum(explained(1:3));