% Running this program will gain the data shown in Figure 5. The running time is estimated to be 3 hours.
clc
clear all
close all
lx=42;%Length of the measurement area x
ly=42;%Length of the measurement area y
down_lim=-lx/2;
up_lim=lx/2;
itv=0.6;%Measurement point interval


%% Import ideal data
A1=load('fig5_new_ideal_saddle.mat');
A2=load('fig5_new_ideal_sink.mat');
A3=load('fig5_new_ideal_source.mat');
A4=load('fig5_new_ideal_NSS_snr.mat');

ideal_saddle_XY=A1.ideal_saddle_XY;
ideal_sink_XY=A2.ideal_sink_XY  ;
ideal_source_XY=A3.ideal_source_XY;
NSS_snr=A4.NSS_snr;

%% Import data with noise


A=load('fig5_new_bxx_600.mat');B=load('fig5_new_bxy_600.mat');C=load('fig5_new_bxz_600.mat');
BXX_noise_total=A.BXX_noise_total;       BXY_noise_total=B.BXY_noise_total;       BXZ_noise_total=C.BXZ_noise_total;
D=load('fig5_new_byy_600.mat');E=load('fig5_new_byz_600.mat');F=load('fig5_new_bzz_600.mat');
BYY_noise_total=D.BYY_noise_total;       BYZ_noise_total=E.BYZ_noise_total;       BZZ_noise_total=F.BZZ_noise_total;
%% Vertical base topology


sum_base=false(1, 15);
sum_rob=false(1, 15);

total_num=600;
for num_flag=1:total_num
    
    for snr_flag=1:15
        SNR_limit=0.8+snr_flag*0.2;%Signal-to-noise ratio limit
        snr_x(snr_flag)=SNR_limit;
        
        
        %% base topology
        
        virtual_M=[0 0 1];
        for i=1:(ly/itv+1)
            for j=1:(lx/itv+1)
                X(i,j)=(j-1)*itv+(-lx/2);
                Y(i,j)=(i-1)*itv+(-ly/2);
                G=[BXX_noise_total(i,j,num_flag) BXY_noise_total(i,j,num_flag) BXZ_noise_total(i,j,num_flag) ;
                    BXY_noise_total(i,j,num_flag) BYY_noise_total(i,j,num_flag) BYZ_noise_total(i,j,num_flag) ;
                    BXZ_noise_total(i,j,num_flag) BYZ_noise_total(i,j,num_flag) BZZ_noise_total(i,j,num_flag) ];
                f=virtual_M*G;
                fx(i,j)=f(1);
                fy(i,j)=f(2);
                norm_fx(i,j)=f(1)/sqrt(f(1)^2+f(2)^2);
                norm_fy(i,j)=f(2)/sqrt(f(1)^2+f(2)^2);
                if NSS_snr(i,j)<SNR_limit
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
        saddle   = collect_points(result.type == 2, X,Y,result);
        
        
        %% Clustering of the identified critical points
        cluster_radius=itv*5;%Clustering radius
        keep_source_XY=[];    keep_sink_XY=[];    keep_saddle_XY =[];
        keep_source_XY = cluster_points_keep_one([source.x source.y],cluster_radius);
        keep_sink_XY = cluster_points_keep_one([sink.x sink.y],cluster_radius);
        keep_saddle_XY = cluster_points_keep_one([saddle.x,saddle.y],cluster_radius*1.3);
        
        base_source_XY=[];    base_sink_XY=[];    base_saddle_XY =[];
        base_source_XY=keep_source_XY;
        base_sink_XY=keep_sink_XY;
        base_saddle_XY=keep_saddle_XY;
        
        
        
        %% Determining the critical point for identifying errors
        alpha=itv*5;%The radius of discrimination
        nearestDist_source=[];  nearestDist_sink=[];    nearestDist_saddle=[];
        idxOutside_source=[];   idxOutside_sink=[];     idxOutside_saddle=[];
        
        [~, nearestDist_source] = knnsearch(ideal_source_XY, base_source_XY);
        idxOutside_source = nearestDist_source > alpha;
        hasOutsidePoint_source = any(idxOutside_source);
        
        
        [~, nearestDist_sink] = knnsearch(ideal_sink_XY, base_sink_XY);
        idxOutside_sink = nearestDist_sink > alpha;
        hasOutsidePoint_sink = any(idxOutside_sink);
        
        
        [~, nearestDist_saddle] = knnsearch(ideal_saddle_XY, base_saddle_XY);
        idxOutside_saddle = nearestDist_saddle > alpha;
        hasOutsidePoint_saddle = any(idxOutside_saddle);
        
        
        %The saddle point is unstable and is not used for the determination of the existence of magnetic sources. Therefore, it is not used for judgment.
        hasOutsidePoint=hasOutsidePoint_source | hasOutsidePoint_sink ;
        single_hasOutsidePoint(snr_flag)=hasOutsidePoint;
        
        %% Vertical robustness topology
        sum_source_XY=[];
        sum_sink_XY=[];
        sum_saddle_XY=[];
        for fl=1:12
            Dh=(fl-1)*30;Ih=80;
            l=cosd(Dh)*cosd(Ih);
            m=sind(Dh)*cosd(Ih);
            n=sind(Ih);
            virtual_M=[l m n];
            for i=1:(ly/itv+1)
                for j=1:(lx/itv+1)
                    
                    G=[BXX_noise_total(i,j,num_flag) BXY_noise_total(i,j,num_flag) BXZ_noise_total(i,j,num_flag) ;
                        BXY_noise_total(i,j,num_flag) BYY_noise_total(i,j,num_flag) BYZ_noise_total(i,j,num_flag) ;
                        BXZ_noise_total(i,j,num_flag) BYZ_noise_total(i,j,num_flag) BZZ_noise_total(i,j,num_flag) ];
                    
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
            keep_source_XY=[];    keep_sink_XY=[];    keep_saddle_XY =[];
            keep_source_XY = cluster_points_keep_one([source.x source.y],cluster_radius);
            keep_sink_XY   = cluster_points_keep_one([sink.x sink.y],cluster_radius);
            keep_saddle_XY = cluster_points_keep_one([saddle.x saddle.y],cluster_radius*1.3);
            
            sum_source_XY=[sum_source_XY;keep_source_XY];
            sum_sink_XY=[sum_sink_XY;keep_sink_XY];
            sum_saddle_XY=[sum_saddle_XY;keep_saddle_XY];
            
        end
     
        
%% Determine how many of the critical points calculated by the rotation of Du are
%% located near the critical points calculated by each vector base topology.

        Discrimination_radius=itv*5;%Search radius
        idx_source=[];count_source=[];idx_source_final=[];source_final=[];
        idx_source = rangesearch(sum_source_XY, base_source_XY, Discrimination_radius);
        count_source = cellfun(@numel, idx_source);
        idx_source_final = find(count_source > 10);
        source_final = base_source_XY(idx_source_final, :);

        idx_sink=[];count_sink=[];idx_sink_final=[];sink_final=[];
        idx_sink = rangesearch(sum_sink_XY, base_sink_XY, Discrimination_radius);
        count_sink = cellfun(@numel, idx_sink);
        idx_sink_final = find(count_sink > 10);
        sink_final = base_sink_XY(idx_sink_final, :);
  
        idx_saddle=[];count_saddle=[];idx_saddle_final=[];saddle_final=[];
        idx_saddle = rangesearch(sum_saddle_XY, base_saddle_XY, Discrimination_radius);
        count_saddle = cellfun(@numel, idx_saddle);
        idx_saddle_final = find(count_saddle > 10);
        saddle_final = base_saddle_XY(idx_saddle_final, :);

        
%% Determine whether there are any errors in the critical points of the retained Vertical robustness topology
        alpha=itv*5;%The radius of discrimination
        nearestDist_source_final=[];    nearestDist_sink_final=[];      nearestDist_saddle_final=[];
        idxOutside_source_final=[];       idxOutside_sink_final=[];     idxOutside_saddle_final=[];
        [~, nearestDist_source_final] = knnsearch( ideal_source_XY,source_final);
        idxOutside_source_final = nearestDist_source_final > alpha;
        hasOutsidePoint_source_final = any(idxOutside_source_final);
        
        [~, nearestDist_sink_final] = knnsearch( ideal_sink_XY,sink_final);
        idxOutside_sink_final = nearestDist_sink_final > alpha;
        hasOutsidePoint_sink_final = any(idxOutside_sink_final);
        
        
        [~, nearestDist_saddle_final] = knnsearch( ideal_saddle_XY,saddle_final);
        idxOutside_saddle_final = nearestDist_saddle_final > alpha;
        hasOutsidePoint_saddle_final = any(idxOutside_saddle_final);

        hasOutsidePoint_final=hasOutsidePoint_source_final | hasOutsidePoint_sink_final ;
        single_hasOutsidePoint_final(snr_flag)=hasOutsidePoint_final;
        
    end
    sum_base=sum_base+single_hasOutsidePoint;
    sum_rob=sum_rob+single_hasOutsidePoint_final;
end
plot(snr_x,sum_base/total_num)
hold on
plot(snr_x,sum_rob/total_num)
legend('base',' rob')
