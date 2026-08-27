clc
clear all
close all
%% import data
lx=40;%Length of the measurement area x
ly=40;%Length of the measurement area y
down_lim=-lx/2;
up_lim=lx/2;
itv=0.8;%Measurement point interval
A=load('fig10_bxx_noise.mat');B=load('fig10_bxy_noise.mat');C=load('fig10_bxz_noise.mat');
bxx_noise=A.bxx_noise;       bxy_noise=B.bxy_noise;       bxz_noise=C.bxz_noise;
D=load('fig10_byy_noise.mat');E=load('fig10_byz_noise.mat');F=load('fig10_bzz_noise.mat');
byy_noise=D.byy_noise;       byz_noise=E.byz_noise;       bzz_noise=F.bzz_noise;
H=load('fig10_mask.mat');
mask=H.mask;
figure

t = tiledlayout(2, 5, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');
for fl=1:10
    %     Change the direction of the virtual magnetic moment and set the virtual magnetic inclination angle to 0°
    Dh=(fl-1)*20+0;Ih=0;
    l=cosd(Dh)*cosd(Ih);
    m=sind(Dh)*cosd(Ih);
    n=sind(Ih);
    virtual_M=[l m n];
    
    for i=1:(ly/itv+1)
        for j=1:(lx/itv+1)
            
            my=-ly/2+(i-1)*itv;mx=-lx/2+(j-1)*itv;X(i,j)=mx;Y(i,j)=my;
            G=[bxx_noise(i,j) bxy_noise(i,j) bxz_noise(i,j) ;
                bxy_noise(i,j) byy_noise(i,j) byz_noise(i,j) ;
                bxz_noise(i,j) byz_noise(i,j) bzz_noise(i,j) ];
            f=virtual_M*G;
            fx(i,j)=f(1);
            fy(i,j)=f(2);
            norm_fx(i,j)=f(1)/sqrt(f(1)^2+f(2)^2);
            norm_fy(i,j)=f(2)/sqrt(f(1)^2+f(2)^2);
            
        end
    end
    
    %     Identify critical points and their types
    rowRadius=1;
    colRadius=1;
    result = struct();    source = struct();    sink   = struct();    saddle = struct();
    result = classify_critical_points_rectangular_ring_with_focus_candidates(X,Y,norm_fx ,norm_fy ,rowRadius,colRadius);
    source = collect_points(result.type == 1,  X,Y,result);
    sink   = collect_points(result.type == -1, X,Y,result);
    saddle   = collect_points(result.type == 2, X,Y,result);
    
    cluster_radius=itv*3;
    keep_source_XY=[];    keep_sink_XY=[];    keep_saddle_XY =[];
    keep_source_XY = cluster_points_keep_one([source.x source.y],cluster_radius);
    keep_sink_XY = cluster_points_keep_one([sink.x sink.y],cluster_radius);
    keep_saddle_XY = cluster_points_keep_one([saddle.x,saddle.y],cluster_radius);
    
    %     Remove the critical points caused by excessive noise.
    k1_flag=1;final_source_XY=[];
    for k1=1:length(keep_source_XY)
        lie=(keep_source_XY(k1,1)-(-lx/2))/itv+1;
        hang=(keep_source_XY(k1,2)-(-ly/2))/itv+1;
        lie=floor(lie);
        hang=floor(hang);
        if (mask(hang,lie)==1)
            final_source_XY(k1_flag,1)=keep_source_XY(k1,1);
            final_source_XY(k1_flag,2)=keep_source_XY(k1,2);
            k1_flag= k1_flag+1;
        end
    end
    
    k2_flag=1;final_sink_XY=[];
    for k2=1:length(keep_sink_XY)
        lie=(keep_sink_XY(k2,1)-(-lx/2))/itv+1;
        hang=(keep_sink_XY(k2,2)-(-ly/2))/itv+1;
        lie=floor(lie);
        hang=floor(hang);
        if (mask(hang,lie)==1)
            final_sink_XY(k2_flag,1)=keep_sink_XY(k2,1);
            final_sink_XY(k2_flag,2)=keep_sink_XY(k2,2);
            k2_flag= k2_flag+1;
        end
    end
    
    k3_flag=1;final_saddle_XY=[];
    for k3=1:length(keep_saddle_XY)
        lie=(keep_saddle_XY(k3,1)-(-lx/2))/itv+1;
        hang=(keep_saddle_XY(k3,2)-(-ly/2))/itv+1;
        lie=floor(lie);
        hang=floor(hang);
        if (mask(hang,lie)==1)
            final_saddle_XY(k3_flag,1)=keep_saddle_XY(k3,1);
            final_saddle_XY(k3_flag,2)=keep_saddle_XY(k3,2);
            k3_flag= k3_flag+1;
        end
    end
    
    
    
    
    
    nexttile
    quiver(X,Y,norm_fx,norm_fy,'k')
    hold on
    scatter(final_source_XY(:,1),final_source_XY(:,2),100,'b','filled');
    hold on
    scatter(final_sink_XY(:,1),final_sink_XY(:,2),100,'s','m','filled');
    hold on
    scatter(final_saddle_XY(:,1),final_saddle_XY(:,2),200,'p','r','filled');
    
    
    xlim([-lx/2 lx/2])
    ylim([-lx/2 lx/2])
    
    
    xlabel('{\it{x}}(m)')
    ylabel('{\it{y}}(m)')
    title([' Du= ' num2str(Dh)  '° '  ' Iu= ' num2str(Ih)  '° '], 'FontWeight','bold')
end

