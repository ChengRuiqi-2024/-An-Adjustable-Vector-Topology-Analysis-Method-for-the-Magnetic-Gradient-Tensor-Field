clear all
clc
close all
natural_linear_file = 'geyer_natural_linear_20m.mat';
rbf_file = 'geyer_rbf_20m.mat';

A = load(natural_linear_file);
B = load(rbf_file);

X = A.X;
Y = A.Y;
valid_mask = A.valid_mask;
FTG_natural = A.FTG_natural;
FTG_linear = A.FTG_linear;
FTG_rbf = B.FTG_rbf;

% RBF: six separate 2-D matrices.
Hxx_rbf = FTG_rbf(:, :, 1);
Hxy_rbf = FTG_rbf(:, :, 2);
Hxz_rbf = FTG_rbf(:, :, 3);
Hyy_rbf = FTG_rbf(:, :, 4);
Hyz_rbf = FTG_rbf(:, :, 5);
Hzz_rbf = FTG_rbf(:, :, 6);


[hang,lie]=size(Hxx_rbf);
NSS = nan(hang, lie);  % Set all to NaN first
valid = isfinite(Hxx_rbf) & isfinite(Hxy_rbf) & isfinite(Hxz_rbf) & ...
    isfinite(Hyy_rbf) & isfinite(Hyz_rbf) & isfinite(Hzz_rbf);

Hxx_rbf(~valid_mask) = 0;
Hxy_rbf(~valid_mask) = 0;
Hxz_rbf(~valid_mask) = 0;
Hyy_rbf(~valid_mask) = 0;
Hyz_rbf(~valid_mask) = 0;
Hzz_rbf(~valid_mask) = 0;

%% Capture window
x_min=min(X(:));
y_min=min(Y(:));

x_left=-2680;x_right=-2260;
y_down=-1340;y_up=-600;

lie_str=(x_left-x_min)/20+1;
lie_end=(x_right-x_min)/20+1;
hang_str=(y_down-y_min)/20+1;
hang_end=(y_up-y_min)/20+1;
X_cut=X(hang_str:hang_end,lie_str:lie_end);
Y_cut=Y(hang_str:hang_end,lie_str:lie_end);

[hang_cut,lie_cut]=size(X_cut);



Hxx_rbf_cut=Hxx_rbf(hang_str:hang_end,lie_str:lie_end);
Hxy_rbf_cut=Hxy_rbf(hang_str:hang_end,lie_str:lie_end);
Hxz_rbf_cut=Hxz_rbf(hang_str:hang_end,lie_str:lie_end);
Hyy_rbf_cut=Hyy_rbf(hang_str:hang_end,lie_str:lie_end);
Hyz_rbf_cut=Hyz_rbf(hang_str:hang_end,lie_str:lie_end);
Hzz_rbf_cut=Hzz_rbf(hang_str:hang_end,lie_str:lie_end);
%% Vertical base topology

virtual_M=[0 0 1];
for i=1:hang_cut
    for j=1:lie_cut

        BXX=Hxx_rbf_cut(i,j);BXY=Hxy_rbf_cut(i,j);BXZ=Hxz_rbf_cut(i,j);BYY=Hyy_rbf_cut(i,j);BYZ=Hyz_rbf_cut(i,j);BZZ=Hzz_rbf_cut(i,j);
        G=[BXX BXY BXZ ;
            BXY BYY BYZ ;
            BXZ BYZ BZZ ];
        f=virtual_M*G;
        fx(i,j)=f(1);
        fy(i,j)=f(2);
        moF(i,j)=f(1)^2+f(2)^2;
        norm_fx_cut(i,j)=f(1)/sqrt(f(1)^2+f(2)^2);
        norm_fy_cut(i,j)=f(2)/sqrt(f(1)^2+f(2)^2);
    end
end

% Identify the critical point
rowRadius=1;
colRadius=1;
result = classify_critical_points_rectangular_ring_with_focus_candidates(X_cut,Y_cut,norm_fx_cut,norm_fy_cut,rowRadius,colRadius);
source = collect_points(result.type == 1,  X_cut,Y_cut,result);
sink   = collect_points(result.type == -1, X_cut,Y_cut,result);
saddle   = collect_points(result.type == 2, X_cut,Y_cut,result);

%cluster
keep_source_XY = cluster_points_keep_one([source.x source.y],20);
keep_sink_XY = cluster_points_keep_one([sink.x sink.y],20);
keep_saddle_XY = cluster_points_keep_one([saddle.x,saddle.y],20);

base_source_XY=keep_source_XY;
base_sink_XY=keep_sink_XY;
base_saddle_XY=keep_saddle_XY;
base_norm_fx_cut=norm_fx_cut;
base_norm_fy_cut=norm_fy_cut;


figure
quiver(X_cut,Y_cut,norm_fx_cut,norm_fy_cut,'k' )
hold on
h3=scatter(base_source_XY(:,1),base_source_XY(:,2),100,'b','filled')
hold on
h4=scatter(base_sink_XY(:,1),base_sink_XY(:,2),100,'s','m','filled')
hold on
h5=scatter(base_saddle_XY(:,1),base_saddle_XY(:,2),200,'p','r','filled')
hold on
h6=scatter(-2560, -1060, 300, [1 0.5 0], 'filled', ...
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.6);
p_in=size(base_source_XY)+size(base_sink_XY)-size(base_saddle_XY);
p_index=p_in(1);
    xlim([x_left x_right])
    ylim([y_down y_up])
  legend([h3,h4,h5,h6],'Source','Sink','Saddle','NSS anomaly', 'FontSize', 12)
title('Vertical base topology of Geyer', 'FontWeight', 'bold', 'FontSize', 16)
xlabel('\itx \rm(m)', 'FontSize', 14);
ylabel('\ity \rm(m)', 'FontSize', 14);
box on
pbaspect([2 3 1])
%% Vertical robustness
sum_source_XY=[];
sum_sink_XY=[];
sum_saddle_XY=[];

for fl=1:12
    
    Dh=(fl-1)*30;Ih=80;
    l=cosd(Dh)*cosd(Ih);
    m=sind(Dh)*cosd(Ih);
    n=sind(Ih);
    virtual_M=[l m n];
    for i=1:hang_cut
        for j=1:lie_cut

            BXX=Hxx_rbf_cut(i,j);BXY=Hxy_rbf_cut(i,j);BXZ=Hxz_rbf_cut(i,j);BYY=Hyy_rbf_cut(i,j);BYZ=Hyz_rbf_cut(i,j);BZZ=Hzz_rbf_cut(i,j);
            G=[BXX BXY BXZ ;
                BXY BYY BYZ ;
                BXZ BYZ BZZ ];
            f=virtual_M*G;
            fx(i,j)=f(1);
            fy(i,j)=f(2);
            moF(i,j)=f(1)^2+f(2)^2;
            norm_fx_cut(i,j)=f(1)/sqrt(f(1)^2+f(2)^2);
            norm_fy_cut(i,j)=f(2)/sqrt(f(1)^2+f(2)^2);
        end
    end
    
    
    rowRadius=1;
    colRadius=1;
    
    result = struct();
    source = struct();
    sink   = struct();
    saddle = struct();
    
    result = classify_critical_points_rectangular_ring_with_focus_candidates(X_cut,Y_cut,norm_fx_cut,norm_fy_cut,rowRadius,colRadius);
    source = collect_points(result.type == 1,  X_cut,Y_cut,result);
    sink   = collect_points(result.type == -1, X_cut,Y_cut,result);
    saddle = collect_points(result.type == 2,  X_cut,Y_cut,result);
    
    keep_source_XY=[];
    keep_sink_XY=[];
    keep_saddle_XY =[];
    keep_source_XY = cluster_points_keep_one([source.x source.y],20);
    keep_sink_XY = cluster_points_keep_one([sink.x sink.y],20);
    keep_saddle_XY = cluster_points_keep_one([saddle.x saddle.y],20);
    
    sum_source_XY=[sum_source_XY;keep_source_XY];
    sum_sink_XY=[sum_sink_XY;keep_sink_XY];    
    sum_saddle_XY=[sum_saddle_XY;keep_saddle_XY];    

    
end
% Determine whether the critical point is stable
idx_source = rangesearch(sum_source_XY, base_source_XY, 80);  
count_source = cellfun(@numel, idx_source);
idx_source_final = find(count_source > 10);
source_final = base_source_XY(idx_source_final, :);


idx_sink = rangesearch(sum_sink_XY, base_sink_XY, 80);  
count_sink = cellfun(@numel, idx_sink);
idx_sink_final = find(count_sink > 10);
sink_final = base_sink_XY(idx_sink_final, :);

idx_saddle = rangesearch(sum_saddle_XY, base_saddle_XY, 80);  
count_saddle = cellfun(@numel, idx_saddle);
idx_saddle_final = find(count_saddle > 10);
saddle_final = base_saddle_XY(idx_saddle_final, :);


figure
hold on
scatter(source_final(:,1),source_final(:,2),100,'b','filled')
hold on
scatter(sink_final(:,1),sink_final(:,2),100,'s','m','filled')
hold on
scatter(saddle_final(:,1),saddle_final(:,2),200,'p','r','filled')
hold on
 scatter(-2560, -1060, 300, [1 0.5 0],'^', 'filled', ...
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.6);

legend('Source','Sink','Saddle','NSS anomaly', 'FontSize', 12)
xlim([x_left x_right])
ylim([y_down y_up])
title('Vertical robustness topology of Geyer', 'FontWeight', 'bold', 'FontSize', 16)
xlabel('\itx \rm(m)', 'FontSize', 14);
ylabel('\ity \rm(m)', 'FontSize', 14);
box on
pbaspect([2 3 1])
% The method for calculating the streamline is not the main focus of this article and is related to subsequent research. 
% Therefore, it will not be made public.
