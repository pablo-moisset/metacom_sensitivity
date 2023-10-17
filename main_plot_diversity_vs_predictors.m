%Run regression_foodweb or regression landscape first, to load
%table_0 or table_1

assert(all(table_0.foodweb_temp==table_1.foodweb_temp))
assert(all(table_0.foodweb_diagonal==table_1.foodweb_diagonal))
assert(all(table_0.MR_a==table_1.MR_a))

toplot = table_0.foodweb_temp >= 0 ; %all true

toplot = toplot & table_0.foodweb_temp == 1 ;     %  0    0.2000    0.4000    0.6000    0.8000    1.0000    1.2000
toplot = toplot & table_0.foodweb_diagonal >= -0.34 ; % -1.0000   -0.8667   -0.7333   -0.6000   -0.4667   -0.3333
%toplot = toplot & table_0.foodweb_diagonal > -0.34 ;
%a = 3000; toplot = toplot & table_0.MR_a == a & table_1.p_gamma<=1+0.75 ;
%a = -1; toplot = toplot & table_0.foodweb_diagonal == a & table_1.p_gamma<=1+0.75 ;
response_var = 'p_gamma'
predictor = 'MR_a'

eval ( ['data_0 = table_0.' response_var '(toplot)'] ) ;
eval ( ['x_values = log10(table_0.' predictor '(toplot))'] ) ;
figure;
subplot(2,1,1) ;
scatter(x_values, data_0) ;
title('A=0') ;
xlabel(predictor) ;
ylabel(response_var);


eval ( ['data_1 = table_1.' response_var '(toplot)'] ) ;
subplot(2,1,2) ;
scatter(x_values,data_1) ;
title('A=1') ;
xlabel(predictor) ;
ylabel(response_var);
set(gcf,'position',[100 100 250 500]) ;
