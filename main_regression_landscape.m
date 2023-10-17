% Computes the main results (panel of regression coefficients)
% Each regression is for fixed values of landscape size and excess factor
% This is a script that requires no parameters, everything is hardcoded
% Assumes the files resulting from post processing are in MOD_landscape/


% These parameters are for axis' scales
L_b_beta.minval  = -0.8; L_b_beta.maxval = 1; 
L_b_gamma.minval = -300; L_b_gamma.maxval = 50; 
L_p_alpha.minval = -0.35; L_p_alpha.maxval = 0.05; 
L_p_beta.minval  = -13; L_p_beta.maxval = 16;
L_p_gamma.minval = -0.1; L_p_gamma.maxval = 0.2; 
L_shannon_alpha.minval = -0.2; L_shannon_alpha.maxval = 0.05; 
L_shannon_beta.minval  = -0.1; L_shannon_beta.maxval = 0.25; 
L_shannon_gamma.minval  = -0.04; L_shannon_gamma.maxval = 0.1; 

% We fix these parameters here
s_richness = 45 ;
landscape_centers = 5 ;
foodweb_connectance = 0.25 ;
landscape_p_exponent = 2 ;

% Create a csv file with experiment data. It is not used here but
% it is a useful format for other analysis
text_output_f = fopen('landscape_panel_data.csv', 'w' ) ;

for landscape_size = 50[10,25,50,100] 
for landscape_excess = 50[5 10 50] 

% Input files for asynchrony=0 and asynchrony=0.5
filename_v0  = sprintf('MOD_landscape/S%d_p4/outputs_p%02dle%02dS%dc%dx%du_var0_*px%d_CON%g',...
                       s_richness, landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent,foodweb_connectance) ;
filename_v05 = sprintf('MOD_landscape/S%d_p4/outputs_p%02dle%02dS%dc%dx%du_var0.5_*px%d_CON%g',...
                       s_richness, landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent,foodweb_connectance) ;

table_0 = load_pp_files(filename_v0);
table_1 = load_pp_files(filename_v05);

response_vars = {'p_alpha',
'p_beta',
'p_gamma' , %minval = -0.08; maxval = 0.045 ;
'shannon_alpha', %minval = -2; maxval = 5.5 ; 
'shannon_beta', %minval = -8.5; maxval = 4 ;
'shannon_gamma', %minval = -3.2; maxval = 2.2 ;
'b_beta' , %minval = -0.5; maxval = 0.2 ;
'b_gamma'} ; %minval = -20; maxval = 80 ;

center=100
for i = 1:numel(response_vars) %for all response variables do a regression and a plot
response_var = response_vars{i}
eval( ['d_response = table_0.' response_var ';'] ) ;
%compute main predictors
predictors = [table_1.foodweb_temp, -table_1.foodweb_diagonal, log10(table_1.MR_a)] ;
predictors = zscore(predictors) ;

X1 = predictors(:,1) ;
X2 = predictors(:,2) ;
X3 = predictors(:,3) ;

%main plus derived predictors
semi_z_predictors = [predictors, X1.*X2, X1.*X3, X2.*X3, X1.*X1, X2.*X2, X3.*X3] ;


% d_response is divided by delta variability, to make it more
% like a derivative
d_response = -2*(eval( ['table_0.' response_var] )- eval( ['table_1.' response_var] )) ;
mdl_d_alt =fitlm(semi_z_predictors, d_response);

%{
%this plots residues, was used for testing only
to_plot = true(size(d_response)) ;
filtered_response = d_response(to_plot) ;
filtered_semi_z_predictors = semi_z_predictors(to_plot,:) ;
mask = [2,8] ; %do not remove these predictors
coeff_to_remove = mdl_d_alt.Coefficients.Estimate(2:end) ;
coeff_to_remove(mask) = 0 ;
partial_residues = filtered_response - sum(filtered_semi_z_predictors.*repmat(coeff_to_remove',sum(to_plot),1),2);
figure
scatter(filtered_semi_z_predictors(:,min(mask))+rand(size(partial_residues,1),1)*0.3,partial_residues) ;
mdl_residue =fitlm(filtered_semi_z_predictors(:,min(mask)),partial_residues);
legend(sprintf("%s",mdl_residue.Coefficients.Estimate(2))) ;
%}

figure
bh=bar(mdl_d_alt.Coefficients.Estimate(2:end), 'FaceColor',[100,170,166]/255) ;
bh.EdgeColor = 'none' ;
ylim([eval(['L_' response_var '.minval']) eval(['L_' response_var '.maxval'])]) ;
hold on ;
h=errorbar(1:size(semi_z_predictors,2), mdl_d_alt.Coefficients.Estimate(2:end),...
           mdl_d_alt.Coefficients.SE(2:end),...
           mdl_d_alt.Coefficients.SE(2:end));
h.LineStyle = 'none' ;
h.Color = 'Black' ;
xlabel('Predictor')
ylabel('Effect size')

set(gcf,'position',[center,center,350,300])
center=center+10;
%export_filename = sprintf('%s_MOD_p%02dS%dc%dx%dpx%d_CON%g',...
%                      response_var,...
%                      landscape_size,s_richness,landscape_centers,...
%                      landscape_excess,landscape_p_exponent,foodweb_connectance);
%saveas(gcf, [export_filename '.fig']);
%close all;
%open([export_filename '.fig']);
%set(gcf,'position',[center,center,350,300])
%saveas(gcf, [export_filename '.tif'], 'tiff');
%tmp = gcf ;
%tmp.Position
%export_fig(['E_' export_filename]);

fprintf(text_output_f, '%s,p%02d,x%d,',response_var,landscape_size,landscape_excess);
fprintf(text_output_f, '%d,', mdl_d_alt.Coefficients.Estimate(2:end)) ;
tmp = sprintf('%d,',mdl_d_alt.Coefficients.SE(2:end));
tmp = tmp(1:(end-1)) ;
fprintf(text_output_f, '%s\n', tmp) ;
%close all;
end
%close all;
end
end
fclose(text_output_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads post processing files with the data for the regression
% experiment pattern uses wild cards to select the appropriate files
function table_all_pp = load_pp_files(experiment_pattern)
   file_pattern = sprintf('%s/pp_*.mat',experiment_pattern) ;
   file_list = dir(file_pattern) ;
   file_table = struct2table(file_list) ;
   file_table = sortrows(file_table, {'folder', 'name'}) ;
   file_list = table2struct(file_table) ; % sort entries by folder+name
   
   %initialize the output structure
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
   
   for i = 1:length(file_list) %read each file and append to arrays in result structure
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
