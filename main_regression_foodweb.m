% Computes the main results (panel of regression coefficients)
% Each regression is for fixed values of foodweb richness and connectance
% This is a script that requires no parameters, everything is hardcoded
% Assumes the files resulting from post processing are in MOD_landscape/


landscape_size = 50;
landscape_centers = 5 ;
landscape_excess = 50 ;
landscape_p_exponent = 2 ;

text_output_f = fopen('foodweb_panel_data.csv', 'w' ) ;

L_b_beta.minval  = -0.8; L_b_beta.maxval = 1; 
L_b_gamma.minval = -300; L_b_gamma.maxval = 50; 
L_p_alpha.minval = -0.35; L_p_alpha.maxval = 0.05; 
L_p_beta.minval  = -4; L_p_beta.maxval = 2;
L_p_gamma.minval = -0.1; L_p_gamma.maxval = 0.2; 
L_shannon_alpha.minval = -0.2; L_shannon_alpha.maxval = 0.05; 
L_shannon_beta.minval  = -0.1; L_shannon_beta.maxval = 25; 
L_shannon_gamma.minval  = -0.04; L_shannon_gamma.maxval = 0.1; 


for s_richness = [30,45,60,75] 
for foodweb_connectance = [0.15 0.2 0.25]
%close all

filename_v0  = sprintf('MOD_landscape/S%d_p4/outputs_p%02dle%02dS%dc%dx%du_var0_*px%d_CON%g',...
                       s_richness, landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent,foodweb_connectance) ;
%filename_v0  = sprintf('MOD_regression_all/outputs_p%02dle%02dS%dc%dx%du_var0_*px%d',...
%                      landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent) ;
filename_v05 = sprintf('MOD_landscape/S%d_p4/outputs_p%02dle%02dS%dc%dx%du_var0.5_*px%d_CON%g',...
                       s_richness, landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent,foodweb_connectance) ;
%filename_v05 = sprintf('MOD_regression_all/outputs_p%02dle%02dS%dc%dx%du_var0.5_*px%d',...
%                      landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent) ;

table_0 = load_pp_files(filename_v0);
table_1 = load_pp_files(filename_v05);

response_vars = {
'p_alpha',...
'p_beta',...
'p_gamma' ,... %minval = -0.08; maxval = 0.045 ;
'shannon_alpha',... %minval = -2; maxval = 5.5 ; 
'shannon_beta',... %minval = -8.5; maxval = 4 ;
'shannon_gamma',... %minval = -3.2; maxval = 2.2 ;
'b_beta' ,... %minval = -0.5; maxval = 0.2 ;
'b_gamma'} ; %minval = -20; maxval = 80 ;


for i=1:numel(response_vars)
response_var = response_vars{i}
%eval( ['d_response = table_0.' response_var ';'] ) ;
predictors = [table_1.foodweb_temp, -table_1.foodweb_diagonal, log10(table_1.MR_a)] ;
predictors = zscore(predictors) ;

X1 = predictors(:,1) ;
X2 = predictors(:,2) ;
X3 = predictors(:,3) ;

semi_z_predictors = [predictors, X1.*X2, X1.*X3, X2.*X3, X1.*X1, X2.*X2, X3.*X3] ;

% d_response is divided by delta variability, to make it more
% like a derivative
%d_response = (eval( ['table_0.' response_var] )- eval( ['table_1.' response_var] )) ;
d_response = -2*(eval( ['table_0.' response_var] )- eval( ['table_1.' response_var] )) ;
mdl_d_alt =fitlm(semi_z_predictors, d_response);

%this plots residues, for testing the code only
%to_plot = true(size(d_response)) ;
%filtered_response = d_response(to_plot) ;
%filtered_semi_z_predictors = semi_z_predictors(to_plot,:) ;
%mask = [1] ; %do not remove these predictors
%coeff_to_remove = mdl_d_alt.Coefficients.Estimate(2:end) ;
%coeff_to_remove(mask) = 0 ;
%partial_residues = filtered_response - sum(filtered_semi_z_predictors.*repmat(coeff_to_remove',sum(to_plot),1),2);
%figure
%scatter(filtered_semi_z_predictors(:,3)+rand(size(partial_residues,1),1)*0.3,partial_residues) ;
%mdl_residue =fitlm(filtered_semi_z_predictors(:,3),partial_residues);
%legend(sprintf("%s",mdl_residue.Coefficients.Estimate(2))) ;

figure
bh=bar(mdl_d_alt.Coefficients.Estimate(2:end), 'FaceColor',[100,170,166]/255) ;
bh.EdgeColor = 'none' ;
ylim([eval(['L_' response_var '.minval']) eval(['L_' response_var '.maxval'])]) ;
%yticks([-0.4:0.2:0.4]) ;
hold on ;
h=errorbar(1:size(semi_z_predictors,2), mdl_d_alt.Coefficients.Estimate(2:end),...
           mdl_d_alt.Coefficients.SE(2:end),...
           mdl_d_alt.Coefficients.SE(2:end));
h.LineStyle = 'none' ;
h.Color = 'Black' ;
%xlabel('Predictor')
%ylabel('Effect size')

%plot([0 10],[0 0],'k')
set(gcf,'position',[100,100,700,300])
set(gcf,'position',[100,100,350,300])
filename = sprintf('%s_MOD_p%02dS%dc%dx%dpx%d_CON%g',...
                   response_var,...
                   landscape_size,s_richness,landscape_centers,...
                   landscape_excess,landscape_p_exponent,foodweb_connectance) ;
savefig([filename '.fig']);
saveas(gcf, [filename '.tif'],'tiff');
fprintf(text_output_f, '%s,S%02d,CON%.2f,',response_var,s_richness,foodweb_connectance);
fprintf(text_output_f, '%d,', mdl_d_alt.Coefficients.Estimate(2:end)) ;
tmp = sprintf('%d,',mdl_d_alt.Coefficients.SE(2:end));
tmp = tmp(1:(end-1)) ;
fprintf(text_output_f, '%s\n', tmp) ;
end
end
end

figure; hist(d_response,20)
%close all
fclose(text_output_f);

%
response_var= 'b_gamma';
eval( ['response_V0 = table_0.' response_var ';'] ) ;
eval( ['response_V1 = table_1.' response_var ';'] ) ;
figure
hist(response_V0(table_0.foodweb_temp==0)) ;
title('V=0 T=0')
figure
hist(response_V1(table_1.foodweb_temp==0)) ;
title('V=1 T=0')
figure
hist(response_V0(table_0.foodweb_temp==1)) ;
title('V=0 T=1')
figure
hist(response_V1(table_1.foodweb_temp==1)) ;
title('V=1 T=1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table_all_pp = load_pp_files(experiment_pattern)
   file_pattern = sprintf('%s/pp_*.mat',experiment_pattern) ;
   file_list = dir(file_pattern) ;
   file_table = struct2table(file_list) ;
   file_table = sortrows(file_table, {'folder', 'name'}) ;
   file_list = table2struct(file_table) ; % sort entries by folder+name
   
   table_all_pp.C = [] ;
   table_all_pp.C_measured = [] ;
   table_all_pp.S = [] ;
   table_all_pp.umig_rate = [] ;
   table_all_pp.sigma = [] ;
   table_all_pp.p_alpha = [] ;
   table_all_pp.p_beta = [] ;
   table_all_pp.p_gamma = [] ;
   table_all_pp.b_beta = [] ;
   table_all_pp.b_gamma = [] ;
   table_all_pp.ext_tsh = [] ;
   table_all_pp.variability = [] ;
   table_all_pp.shannon_alpha = [] ;
   table_all_pp.shannon_beta = [] ;
   table_all_pp.shannon_gamma = [] ;
   table_all_pp.wetness = [] ;
   table_all_pp.foodweb_temp = [] ;
   table_all_pp.foodweb_diagonal = [] ;
   table_all_pp.MR_a = [] ;
   
   for i = 1:length(file_list) 
       fprintf('File nbr. %d\n',i)
       fname = [file_list(i).folder '/' file_list(i).name] ;
       load(fname) ;
 %      assert(all(all_p_gamma<=1))
       table_all_pp.C = [table_all_pp.C; all_C] ;
       table_all_pp.C_measured = [table_all_pp.C_measured; all_C_measured] ;
       table_all_pp.S = [table_all_pp.S; all_S] ;
       table_all_pp.umig_rate = [table_all_pp.umig_rate; all_umig_rate] ;
       table_all_pp.sigma = [table_all_pp.sigma; all_sigma] ;
       %   total_biomass = [table_all_pp.] ;
       table_all_pp.p_alpha = [table_all_pp.p_alpha; all_p_alpha] ;
       table_all_pp.p_beta = [table_all_pp.p_beta; all_p_beta] ;
       table_all_pp.p_gamma = [table_all_pp.p_gamma; all_p_gamma] ;
       table_all_pp.b_beta = [table_all_pp.b_beta; all_b_beta] ;
       table_all_pp.b_gamma = [table_all_pp.b_gamma; all_b_gamma] ;
       table_all_pp.ext_tsh = [table_all_pp.ext_tsh; all_ext_tsh] ;
       table_all_pp.variability = [table_all_pp.variability; all_variability] ;
       table_all_pp.shannon_alpha = [table_all_pp.shannon_alpha; all_shannon_alpha] ;
       table_all_pp.shannon_beta = [table_all_pp.shannon_beta; all_shannon_beta] ;
       table_all_pp.shannon_gamma = [table_all_pp.shannon_gamma; all_shannon_gamma] ;
       tmp = [] ;
       for sample=1:numel(M)
          tmp = [tmp; M{1,sample}.len_on_period'] ;
       end
       table_all_pp.wetness = [table_all_pp.wetness; tmp] ;

       table_all_pp.foodweb_temp = [table_all_pp.foodweb_temp;repmat(MASTER_P.temperature,numel(M),1)] ;
       table_all_pp.foodweb_diagonal = [table_all_pp.foodweb_diagonal;repmat(MASTER_P.diagonal,numel(M),1)] ;
       table_all_pp.MR_a = [table_all_pp.MR_a;repmat(METACOM_SIM_P.mig_rate_p(1),numel(M),1)] ;
   end
end
