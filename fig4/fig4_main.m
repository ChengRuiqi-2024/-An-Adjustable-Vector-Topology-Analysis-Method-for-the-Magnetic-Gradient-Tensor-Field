clc
clear all
close all

lx=42;%Length of the measurement area x
ly=42;%Length of the measurement area y
down_lim=-lx/2;
up_lim=lx/2;
itv=0.6;%Measurement point interval


%% Import ideal data
A1=load('fig4_new_ideal_saddle.mat');
A2=load('fig4_new_ideal_sink.mat');
A3=load('fig4_new_ideal_source.mat');
A4=load('fig4_new_ideal_NSS_snr.mat');

ideal_saddle_XY=A1.ideal_saddle_XY;
ideal_sink_XY=A2.ideal_sink_XY  ;
ideal_source_XY=A3.ideal_source_XY;
NSS_snr=A4.NSS_snr;



%% Import data with noise

SNR_limit=1.4;%%Signal-to-noise ratio limit
A=load('fig4_bxx_noise.mat');B=load('fig4_bxy_noise.mat');C=load('fig4_bxz_noise.mat');
bxx_noise=A.bxx_noise;       bxy_noise=B.bxy_noise;       bxz_noise=C.bxz_noise;
D=load('fig4_byy_noise.mat');E=load('fig4_byz_noise.mat');F=load('fig4_bzz_noise.mat');
byy_noise=D.byy_noise;       byz_noise=E.byz_noise;       bzz_noise=F.bzz_noise;
%% Vertical base topology
virtual_M=[0 0 1];
for i=1:(ly/itv+1)
    for j=1:(lx/itv+1)
        X(i,j)=(j-1)*itv+(-lx/2);
        Y(i,j)=(i-1)*itv+(-ly/2);
        
        G=[bxx_noise(i,j) bxy_noise(i,j) bxz_noise(i,j) ;
            bxy_noise(i,j) byy_noise(i,j) byz_noise(i,j) ;
            bxz_noise(i,j) byz_noise(i,j) bzz_noise(i,j) ];
        
        f=virtual_M*G;
        fx(i,j)=f(1);
        fy(i,j)=f(2);
        norm_fx(i,j)=f(1)/sqrt(f(1)^2+f(2)^2);
        norm_fy(i,j)=f(2)/sqrt(f(1)^2+f(2)^2);
        
        
        if NSS_snr(i,j)<SNR_limit% Remove data that does not meet the signal-to-noise ratio criteria
            norm_fx_dis(i,j)=-1;
            norm_fy_dis(i,j)=-1;
            mask(i,j)=0;
        else
            norm_fx_dis(i,j)=norm_fx(i,j);
            norm_fy_dis(i,j)=norm_fy(i,j);
            mask(i,j)=1;
        end
    end
end
rowRadius=1;
colRadius=1;
result = classify_critical_points_rectangular_ring_with_focus_candidates(X,Y,norm_fx_dis,norm_fy_dis,rowRadius,colRadius);
source = collect_points(result.type == 1,  X,Y,result);
sink   = collect_points(result.type == -1, X,Y,result);
saddle   = collect_points(result.type == 2, X,Y,result);


%% Clustering of the identified critical points
cluster_radius=itv*5;%Clustering radius
keep_source_XY = cluster_points_keep_one([source.x source.y],cluster_radius);
keep_sink_XY = cluster_points_keep_one([sink.x sink.y],cluster_radius);
keep_saddle_XY = cluster_points_keep_one([saddle.x,saddle.y],cluster_radius*1.3);



base_source_XY=keep_source_XY;
base_sink_XY=keep_sink_XY;
base_saddle_XY=keep_saddle_XY;


figure
quiver(X,Y,norm_fx_dis,norm_fy_dis)
hold on
f3_1=scatter(base_source_XY(:,1),base_source_XY(:,2),100,'b','filled');
hold on
f3_2=scatter(base_sink_XY(:,1),base_sink_XY(:,2),100,'s','m','filled');
hold on
f3_3=scatter(base_saddle_XY(:,1),base_saddle_XY(:,2),200,'p','r','filled');



%% Determining the critical point for identifying errors
alpha=itv*5;%The radius of discrimination
[~, nearestDist_source] = knnsearch(ideal_source_XY, base_source_XY);
idxOutside_source = nearestDist_source > alpha;
base_error_source_XY = base_source_XY(idxOutside_source, :);
hasOutsidePoint_source = any(idxOutside_source);


[~, nearestDist_sink] = knnsearch(ideal_sink_XY, base_sink_XY);
idxOutside_sink = nearestDist_sink > alpha;
base_error_sink_XY = base_sink_XY(idxOutside_sink, :);
hasOutsidePoint_sink = any(idxOutside_sink);


[~, nearestDist_saddle] = knnsearch(ideal_saddle_XY, base_saddle_XY);
idxOutside_saddle = nearestDist_saddle > alpha;
base_error_saddle_XY = base_saddle_XY(idxOutside_saddle, :);
hasOutsidePoint_saddle = any(idxOutside_saddle);

hasOutsidePoint=hasOutsidePoint_source | hasOutsidePoint_sink | hasOutsidePoint_saddle;
if hasOutsidePoint == 1
    disp('Vertical base topology has error points')
else
    disp('Vertical base topology has no error critical point')
end
%
% f3_4=scatter(base_error_source_XY(:,1),base_error_source_XY(:,2),100,'k','filled');
% hold on
% f3_5=scatter(base_error_sink_XY(:,1),base_error_sink_XY(:,2),100,'s','k','filled');
% hold on
% f3_6=scatter(base_error_saddle_XY(:,1),base_error_saddle_XY(:,2),200,'p','k','filled');

% legend([f3_1,f3_2,f3_3,f3_4,f3_5,f3_6],'source','sink','saddle','error source','error sink','error saddle')
legend([f3_1,f3_2,f3_3],'source','sink','saddle')
title('The calculation results of the vertical base topology show that the error points cannot be eliminated.')
xlim([-21 21])
ylim([-21 21])
%% Vertical robustness topology
sum_source_XY=[];
sum_sink_XY=[];
sum_saddle_XY=[];

for fl=1:12
    %Change the direction of the virtual magnetic moment
    Dh=(fl-1)*30;Ih=80;
    l=cosd(Dh)*cosd(Ih);
    m=sind(Dh)*cosd(Ih);
    n=sind(Ih);
    virtual_M=[l m n];
    for i=1:(ly/itv+1)
        for j=1:(lx/itv+1)
            
            G=[bxx_noise(i,j) bxy_noise(i,j) bxz_noise(i,j) ;
                bxy_noise(i,j) byy_noise(i,j) byz_noise(i,j) ;
                bxz_noise(i,j) byz_noise(i,j) bzz_noise(i,j) ];
            
            f=virtual_M*G;
            fx(i,j)=f(1);
            fy(i,j)=f(2);
            norm_fx(i,j)=f(1)/sqrt(f(1)^2+f(2)^2);
            norm_fy(i,j)=f(2)/sqrt(f(1)^2+f(2)^2);
            
            if NSS_snr(i,j)<SNR_limit% Remove data that does not meet the signal-to-noise ratio criteria
                norm_fx_dis(i,j)=-1;
                norm_fy_dis(i,j)=-1;
            else
                norm_fx_dis(i,j)=norm_fx(i,j);
                norm_fy_dis(i,j)=norm_fy(i,j);
            end
            
        end
    end
    
    
    rowRadius=1;
    colRadius=1;
    
    result = struct();    source = struct();    sink   = struct();    saddle = struct();
    result = classify_critical_points_rectangular_ring_with_focus_candidates(X,Y,norm_fx_dis,norm_fy_dis,rowRadius,colRadius);
    source = collect_points(result.type == 1,  X,Y,result);
    sink   = collect_points(result.type == -1, X,Y,result);
    saddle = collect_points(result.type == 2,  X,Y,result);
    
    %% Clustering of the identified critical points
    cluster_radius=itv*5;%Clustering radius
    %Sum up the results of each calculation
    keep_source_XY=[];    keep_sink_XY=[];    keep_saddle_XY =[];
    keep_source_XY = cluster_points_keep_one([source.x source.y],cluster_radius);
    keep_sink_XY   = cluster_points_keep_one([sink.x sink.y],cluster_radius);
    keep_saddle_XY = cluster_points_keep_one([saddle.x saddle.y],cluster_radius*1.3);
    
    sum_source_XY=[sum_source_XY;keep_source_XY];
    sum_sink_XY=[sum_sink_XY;keep_sink_XY];
    sum_saddle_XY=[sum_saddle_XY;keep_saddle_XY];
    
    
end

figure

s1=scatter(sum_source_XY(:,1),sum_source_XY(:,2),50,'k', 'LineWidth', 1.5 );
hold on
s2=scatter(sum_sink_XY(:,1),sum_sink_XY(:,2),50,'s','k', 'LineWidth', 1.5 );
hold on
s3=scatter(sum_saddle_XY(:,1),sum_saddle_XY(:,2),100,'p','k', 'LineWidth', 1.5 );
xlim([-21 21])
ylim([-21 21])


%% Determine how many of the critical points calculated by the rotation of Du are
%% located near the critical points calculated by each vector base topology.

Discrimination_radius=itv*5;%Search radius
idx_source = rangesearch(sum_source_XY, base_source_XY, Discrimination_radius);
count_source = cellfun(@numel, idx_source);
idx_source_final = find(count_source > 10);%If the number is greater than 8, then the corresponding vector base topology critical point will be retained.
source_final = base_source_XY(idx_source_final, :);
idx_remain_source_final = [idx_source{idx_source_final}];
remain_source_final = sum_source_XY(idx_remain_source_final, :);
hold on
s4=scatter(remain_source_final(:,1),remain_source_final(:,2),100,'b','filled');

idx_remove_source_final = find(count_source <= 10);
source_remove = base_source_XY(idx_remove_source_final, :);



idx_sink = rangesearch(sum_sink_XY, base_sink_XY, Discrimination_radius);
count_sink = cellfun(@numel, idx_sink);
idx_sink_final = find(count_sink > 10);
sink_final = base_sink_XY(idx_sink_final, :);
idx_remain_sink_final = [idx_sink{idx_sink_final}];
remain_sink_final = sum_sink_XY(idx_remain_sink_final, :);
hold on
s5=scatter(remain_sink_final(:,1),remain_sink_final(:,2),100,'s','m','filled');
idx_remove_sink_final = find(count_sink <= 10);
sink_remove = base_sink_XY(idx_remove_sink_final, :);


Discrimination_radius=Discrimination_radius*1.5;
idx_saddle = rangesearch(sum_saddle_XY, base_saddle_XY, Discrimination_radius);
count_saddle = cellfun(@numel, idx_saddle);
idx_saddle_final = find(count_saddle > 10);
saddle_final = base_saddle_XY(idx_saddle_final, :);
idx_remain_saddle_final = [idx_saddle{idx_saddle_final}];
remain_saddle_final = sum_saddle_XY(idx_remain_saddle_final, :);
hold on
s6=scatter(remain_saddle_final(:,1),remain_saddle_final(:,2),200,'p','r','filled');
idx_remove_saddle_final = find(count_saddle <= 10);
saddle_remove = base_saddle_XY(idx_remove_saddle_final, :);

legend([s4,s5,s6,s1,s2,s3],'Retained source cluster','Retained sink cluster','Retained saddle cluster','Excluded source cluster','Excluded sink cluster','Excluded saddle cluster','FontSize',12)
% legend([s1,s2,s3,s4,s5,s6],'Excluded source cluster','Excluded sink cluster','Excluded saddle cluster','Retained source cluster','Retained sink cluster','Retained saddle cluster','FontSize',12)
title('Accumulated critical points','FontSize',18,'FontWeight','bold')
xlabel('{\it{x}}(m)','FontSize',14)
ylabel('{\it{y}}(m)','FontSize',14)
xlim([-21 21])
ylim([-21 21])
box on
grid on


figure
% quiver(X ,Y ,base_norm_fx_dis ,base_norm_fy_dis  )
hold on
fin5=scatter(source_remove(:,1),source_remove(:,2),50,'k', 'LineWidth', 1.5 );
hold on
fin6=scatter(sink_remove(:,1),sink_remove(:,2),50,'s','k' , 'LineWidth', 1.5);
hold on
fin7=scatter(saddle_remove(:,1),saddle_remove(:,2),100,'p','k' , 'LineWidth', 1.5);
hold on
fin2=scatter(source_final(:,1),source_final(:,2),100,'b','filled');
hold on
fin3=scatter(sink_final(:,1),sink_final(:,2),100,'s','m','filled');
hold on
fin4=scatter(saddle_final(:,1),saddle_final(:,2),200,'p','r','filled');
legend([fin2,fin3,fin4,fin5,fin6,fin7],'Source','Sink','Saddle','Excluded source','Excluded sink','Excluded saddle','FontSize',12)
title('Vertical robustness topology','FontSize',18,'FontWeight','bold')
xlabel('{\it{x}}(m)','FontSize',14)
ylabel('{\it{y}}(m)','FontSize',14)
xlim([-21 21])
ylim([-21 21])
box on
grid on

%% Determine whether there are any errors in the critical points of the retained Vertical robustness topology
alpha=itv*5;%The radius of discrimination
[~, nearestDist_source_final] = knnsearch( ideal_source_XY,source_final);
idxOutside_source_final = nearestDist_source_final > alpha;
final_error_source_XY = source_final(idxOutside_source_final, :);
hasOutsidePoint_source_final = any(idxOutside_source_final);

[~, nearestDist_sink_final] = knnsearch( ideal_sink_XY,sink_final);
idxOutside_sink_final = nearestDist_sink_final > alpha;
final_error_sink_XY = sink_final(idxOutside_sink_final, :);
hasOutsidePoint_sink_final = any(idxOutside_sink_final);


[~, nearestDist_saddle_final] = knnsearch( ideal_saddle_XY,saddle_final);
idxOutside_saddle_final = nearestDist_saddle_final > alpha;
final_error_saddle_XY = saddle_final(idxOutside_saddle_final, :);
hasOutsidePoint_saddle_final = any(idxOutside_saddle_final);

hasOutsidePoint_final=hasOutsidePoint_source_final | hasOutsidePoint_sink_final | hasOutsidePoint_saddle_final;
if hasOutsidePoint_final == 1
    disp('Vertical robustness topology has error points')
else
    disp('Vertical robustness topology has no error critical point')
end

