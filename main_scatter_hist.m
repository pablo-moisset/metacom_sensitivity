% main script to draw histograms of a response variable after filtering
% based on predictor values. table_0 and table_1 must be defined as in
% main_regression_foodweb.m in this workspace. The easiest way to do it is to run
% main_regression_foodweb.m and after that run this script
% This script also plots the scatter plots for response varaibles for
% two values of asynchrony

assert(all(table_0.foodweb_temp==table_1.foodweb_temp))
assert(all(table_0.foodweb_diagonal==table_1.foodweb_diagonal))
assert(all(table_0.MR_a==table_1.MR_a))

toplot = table_0.foodweb_temp >= 0 ; %all true

%toplot = toplot & table_0.foodweb_temp == 1.2 ;
%toplot = toplot & table_0.foodweb_diagonal > -0.34 ;
a = 3000; toplot = toplot & table_0.MR_a == a & table_1.p_gamma<=1+0.75 ;
%a = -1; toplot = toplot & table_0.foodweb_diagonal == a & table_1.p_gamma<=1+0.75 ;
response_var = 'p_gamma'

f_heatmap = openfig('R_two_patch_p_S45_CON0.2.fig') ;
data_heatmap = f_heatmap.Children.Children.CData ;
data_XData = f_heatmap.Children.Children.XData ;
data_YData = f_heatmap.Children.Children.YData ;
colors_heatmap = RGB_list() ;

eval ( ['data_0 = table_0.' response_var '(toplot)'] ) ;

tmp = table_0.foodweb_temp(toplot) ;
ii = lin_index(tmp,min(tmp),max(tmp),size(data_heatmap,1)) ;
tmp = table_0.foodweb_diagonal(toplot) ;
jj = lin_index(tmp,min(tmp),max(tmp),size(data_heatmap,2)) ;

tmp = [] ;
for k = 1:numel(ii)
    tmp = [tmp;data_heatmap(ii(k),jj(k))] ;
end
ii = lin_index(tmp,min(tmp),max(tmp),size(colors_heatmap,1)) ;

eval ( ['data_1 = table_1.' response_var '(toplot)'] ) ;

%scatter(data_0,data_1,7,colors_heatmap(ii,:),'filled') ;

name = sprintf('L%d_a%d_F%d',landscape_size, a, landscape_excess) ;
figure;hist(data_0) ;
figure;hist(data_1) ;
f=figure('name', name);

colors_dots = colors_heatmap(ii,:) ;
p = randperm(numel(data_0)) ;
scatter(data_0(p),data_1(p),10,colors_dots(p),'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4) ;

f.Position =  [100 500 250 200 ] ; %for consistency

r=refline(1,0);
r.LineStyle = '--' ;
r.Color = [0.5 0.5 0.5];

saveas(f,[name '_' response_var '.fig']);
saveas(f,[name '_' response_var '.pdf']);

toconsider = toplot & (table_1.p_gamma >= 0.5) ; 

diagonals = table_1.foodweb_diagonal(toconsider) ;
temps = table_1.foodweb_temp(toconsider) ;

ndiagonals = table_1.foodweb_diagonal(~toconsider) ;
ntemps = table_1.foodweb_temp(~toconsider) ;

