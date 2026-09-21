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

% RBF: six separate 2-D matrices.
Hxx_rbf = FTG_rbf(:, :, 1);
Hxy_rbf = FTG_rbf(:, :, 2);
Hxz_rbf = FTG_rbf(:, :, 3);
Hyy_rbf = FTG_rbf(:, :, 4);
Hyz_rbf = FTG_rbf(:, :, 5);
Hzz_rbf = FTG_rbf(:, :, 6);



[hang,lie]=size(Hxx_rbf);
NSS = nan(hang, lie);  
valid = isfinite(Hxx_rbf) & isfinite(Hxy_rbf) & isfinite(Hxz_rbf) & ...
    isfinite(Hyy_rbf) & isfinite(Hyz_rbf) & isfinite(Hzz_rbf);
for i=1:hang
    for j=1:lie
        if ~valid(i, j)
            continue   
        end
        BXX=Hxx_rbf(i,j);BXY=Hxy_rbf(i,j);BXZ=Hxz_rbf(i,j);BYY=Hyy_rbf(i,j);BYZ=Hyz_rbf(i,j);BZZ=Hzz_rbf(i,j);
        G=[BXX BXY BXZ ;
            BXY BYY BYZ ;
            BXZ BYZ BZZ ];
        
        % Caculate NSS
        [V, D] = eig(G);
        eigenvalues = diag(D);
        [eigenvalues_sorted, index] = sort(eigenvalues, 'descend');
        lam1=eigenvalues_sorted(1);
        lam2=eigenvalues_sorted(2);
        lam3=eigenvalues_sorted(3);
        NSS(i,j)=sqrt( -lam2^2-lam1*lam3 );
        

    end
end
contourf(X, Y, NSS, 10, 'LineColor', 'none');
cb = colorbar;
title(cb, 'nT m^{-1}', 'Interpreter', 'tex', 'FontSize', 14);
%% vertical
x_min=min(X(:));
y_min=min(Y(:));

x_left_vertical=-2680;x_right_vertical=-2260;
y_down_vertical=-1340;y_up_vertical=-600;

lie_str_vertical=(x_left_vertical-x_min)/20+1;
lie_end_vertical=(x_right_vertical-x_min)/20+1;
hang_str_vertical=(y_down_vertical-y_min)/20+1;
hang_end_vertical=(y_up_vertical-y_min)/20+1;
X_vertical_cut=X(hang_str_vertical:hang_end_vertical,lie_str_vertical:lie_end_vertical);
Y_vertical_cut=Y(hang_str_vertical:hang_end_vertical,lie_str_vertical:lie_end_vertical);
mask_vertical = true(size(X_vertical_cut));
mask_vertical(2:end-1, 2:end-1) = false;  
X_side_vertical = X_vertical_cut(mask_vertical);
Y_side_vertical = Y_vertical_cut(mask_vertical);
hold on
s1=scatter(X_side_vertical,Y_side_vertical,8,'r','filled')
% legend( [s1],'area of vertial topology' , 'FontSize', 10)
axis equal
title('NSS of geyer', 'FontWeight', 'bold', 'FontSize', 16);
xlabel('\it{x} \rm(m)', 'FontSize', 14);
ylabel('\it{y} \rm(m)', 'FontSize', 14);

