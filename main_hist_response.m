%main script to draw histograms of a response variable after filtering
%based on predictor values. table_0 and table_1 must be defined as in
%regression_foodweb.m

assert(all(table_0.foodweb_temp==table_1.foodweb_temp))
assert(all(table_0.foodweb_diagonal==table_1.foodweb_diagonal))
assert(all(table_0.MR_a==table_1.MR_a))

toplot = table_0.foodweb_temp >= 0 ; %all true

%toplot = toplot & table_0.foodweb_temp == 1.2 ;
%toplot = toplot & table_0.foodweb_diagonal > -0.34 ;
%a = 30; toplot = toplot & table_0.MR_a == a & table_1.p_gamma<=1+0.75 ;

%a = -1.000000000000000
%a = -0.600000000000000
a = -0.333333000000000
toplot = toplot & table_0.foodweb_diagonal == a & table_1.p_gamma<=1+0.75 ;

response_var = 'p_gamma'

f_heatmap = openfig('R_two_patch_p_S75_CON0.25.fig') ;
data_heatmap = f_heatmap.Children.Children.CData ;
data_XData = f_heatmap.Children.Children.XData ;
data_YData = f_heatmap.Children.Children.YData ;
colors_heatmap = RGB_list() ;

eval ( ['data_0 = table_0.' response_var '(toplot)'] ) ;

tmp = table_0.foodweb_temp(toplot) ;
jj = lin_index(tmp,min(data_XData),max(data_XData),size(data_heatmap,2)) ;
tmp = table_0.foodweb_diagonal(toplot) ;
ii = lin_index(tmp, min(data_YData),max(data_YData),size(data_heatmap,1)) ;

tmp = [] ;
for k = 1:numel(ii)
    tmp = [tmp;data_heatmap(ii(k),jj(k))] ;
end
values_heatmap = tmp ;

eval ( ['data_1 = table_1.' response_var '(toplot)'] ) ;
%data_1 = 2*(table_1.p_gamma(toplot)-table_0.p_gamma(toplot))
%data_1 = table_1.p_alpha(toplot) ;
%data_1 = table_1.p_beta(toplot) ;

%scatter(data_0,data_1,7,colors_heatmap(ii,:),'filled') ;

name = sprintf('L%d_diag%d_F%d',landscape_size, a, landscape_excess) ;
%figure;hist(data_0) ;
%figure;hist(data_1) ;
f=figure('name', name);

color_idx = lin_index(tmp, 0,1 ,size(colors_heatmap,1)) ;

colors_dots = colors_heatmap(color_idx,:) ;
p = randperm(numel(data_0)) ;
%scatter(data_0(p),data_1(p),100,colors_dots(p,:),'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4) ;
scatter(data_0(p),data_1(p),10,colors_dots(p,:),'MarkerFaceAlpha',0,'MarkerEdgeAlpha',.4) ;

%xlabel('$\mathcal{P}_\gamma$ for $V=0$', 'interpreter','latex') ;
%ylabel('$\mathcal{P}_\gamma$ for $V=0.5$', 'interpreter','latex') ;
xlim([0 1]) ;
ylim([0 1]) ;
f.Position =  [100 500 250 200 ] ;

r=refline(1,0);
r.LineStyle = '--' ;
r.Color = [0.5 0.5 0.5];

data_v0 = data_0(p) ;
values_heatmap_p = values_heatmap(p);
%data_v2 = table_0.foodweb_temp(toplot) ;
%data_response = table_1.p_gamma(toplot) ;
%figure; plot3(data_v1,data_v2,data_response) ;

saveas(f,[name '_' response_var '.fig']);
saveas(f,[name '_' response_var '.pdf']);

toconsider = toplot & (table_1.p_gamma >= 0.5) ; 

diagonals = table_1.foodweb_diagonal(toconsider) ;
temps = table_1.foodweb_temp(toconsider) ;

ndiagonals = table_1.foodweb_diagonal(~toconsider) ;
ntemps = table_1.foodweb_temp(~toconsider) ;

