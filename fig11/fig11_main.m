clear all
clc
close all
natural_linear_file = 'geyer_natural_linear_20m.mat';
rbf_file = 'geyer_rbf_20m.mat';
A = load(natural_linear_file);
B = load(rbf_file);

X = A.X;
Y = A.Y;
FTG_rbf = B.FTG_rbf;
valid_mask = A.valid_mask;
% RBF: six separate 2-D matrices.
Hxx_rbf = FTG_rbf(:, :, 1);
Hxy_rbf = FTG_rbf(:, :, 2);
Hxz_rbf = FTG_rbf(:, :, 3);
Hyy_rbf = FTG_rbf(:, :, 4);
Hyz_rbf = FTG_rbf(:, :, 5);
Hzz_rbf = FTG_rbf(:, :, 6);


Hxx_rbf(~valid_mask) = 0.11;
Hxy_rbf(~valid_mask) = 0.12;
Hxz_rbf(~valid_mask) = 0.13;
Hyy_rbf(~valid_mask) = 0.14;
Hyz_rbf(~valid_mask) = 0.15;
Hzz_rbf(~valid_mask) = 0.16;

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


%% Horizontal dynamic topology
figure
t = tiledlayout(4,5, ...
    'TileSpacing','compact', ... 
    'Padding','compact');         


for fl=1:19
    %     Change the direction of the virtual magnetic moment and set the virtual magnetic inclination angle to 0°
    Dh=(fl-1)*10;Ih=0;
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
    result=struct();source =struct();sink =struct();saddle =struct();focus =struct();wn=struct();
    
    result = classify_critical_points_rectangular_ring_with_focus_candidates(X_cut,Y_cut,norm_fx_cut,norm_fy_cut,rowRadius,colRadius);
    source = collect_points(result.type == 1,  X_cut,Y_cut,result);
    sink   = collect_points(result.type == -1, X_cut,Y_cut,result);
    saddle = collect_points(result.type == 2,  X_cut,Y_cut,result);
    focus = collect_points(result.type == 3,  X_cut,Y_cut,result);

    
    keep_source_XY = [];keep_sink_XY = [];keep_saddle_XY = [];keep_focus_XY = [];
    keep_source_XY = cluster_points_keep_one([source.x source.y],20);
    keep_sink_XY = cluster_points_keep_one([sink.x sink.y],20);
    keep_saddle_XY = cluster_points_keep_one([saddle.x saddle.y],20);
    keep_focus_XY = cluster_points_keep_one([focus.x focus.y],20);
    
    count_saddle(fl)=size(keep_saddle_XY,1)
    Dh_line(fl)=Dh;
    nexttile

    quiver(X_cut,Y_cut,norm_fx_cut,norm_fy_cut )
    hold on
    scatter(keep_source_XY(:,1),keep_source_XY(:,2),100,'b','filled')
    hold on
    scatter(keep_sink_XY(:,1),keep_sink_XY(:,2),100,'s','m','filled')
    hold on
    scatter(keep_saddle_XY(:,1),keep_saddle_XY(:,2),200,'p','r','filled')


    
  
    total_s_XY=[];
    total_s_XY=[keep_sink_XY;keep_source_XY];
    
    keep = true(size(keep_focus_XY, 1), 1);  % Each focus point is initially set to remain as it is.
    for flag_f=1:length(keep_focus_XY(:,1))
        for flag_s=1:length(total_s_XY(:,1))
            d=norm(keep_focus_XY(flag_f,:)-total_s_XY(flag_s,:));
            if d<40
            keep(flag_f) = false; %Delete the focal points that are close to the sinks and the sources.
            end
        end
    end
    keep_focus_XY = keep_focus_XY(keep, :);

    
    hold on
    scatter(keep_focus_XY(:,1),keep_focus_XY(:,2),100,'k','filled')
    subtitle([ ' Dh=' num2str(Dh)])
    xlim([x_left x_right])
    ylim([y_down y_up])
    
    
end

count_saddle2=[ count_saddle count_saddle(1)];

for i=1:length(count_saddle)
    p_behavior(i)=count_saddle2(i+1)-count_saddle(i);
end

figure
plot(Dh_line, p_behavior, '-o', ...
    'Color', 'r', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', 'r');

yticks(-2:3)

yticklabels({ ...
    'Two merging pairs', ...      % −2
    'One merging pair', ...       % −1
    'Unchanged', ...         %  0
    'One splitting pair', ...     %  1
    'Two splitting pairs', ...    %  2
    'Three splitting pairs'});    %  3
% 
 title('Critical point bifurcation of Geyer', 'FontWeight', 'bold', 'FontSize', 14);
   xlabel('Virtual magnetic declination (°)')

