response_var = 'p_alpha' ;
%response_var = 'p_gamma' ;
%response_var = 'b_gamma' ;
%response_var = 'b_beta_div' ;
%response_var = 'b_beta_sub' ;
%response_var = 'p_beta_div' ;
%response_var = 'p_beta_sub' ;

landscape_size = 50 ;
landscape_centers = 5 ;
landscape_excess = 50 ;
landscape_p_exponent = 2 ;
s_richness = 60 ;
foodweb_connectance = 0.2 ;

filename_v0  = sprintf('MOD_regression/outputs_p%02dle%02dS%dc%dx%du_var0_*px%d_CON%g',...
                      landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent,foodweb_connectance) ;
filename_v05 = sprintf('MOD_regression/outputs_p%02dle%02dS%dc%dx%du_var0.5_*px%d_CON%g',...
                      landscape_size,2*landscape_size,s_richness,landscape_centers,landscape_excess,landscape_p_exponent,foodweb_connectance) ;

%table_0 = load_pp_files(filename_v0);
%table_1 = load_pp_files(filename_v05);

d_response = eval( ['table_0.' response_var ';'] ) ;
%d_response = eval( ['table_1.' response_var ';'] ) ;
%d_response = eval( ['table_0.' response_var] ) - eval( ['table_1.' response_var] ) ;
%d_response = (eval( ['table_0.' response_var] ) - eval( ['table_1.' response_var] ))./eval( ['table_0.' response_var] ) ;

predictors = [table_1.foodweb_temp, -table_1.foodweb_diagonal, log10(table_1.MR_a)] ;
X1 = predictors(:,1) ;
X2 = predictors(:,2) ;
X3 = predictors(:,3) ;
predictors = [predictors, X1.*X2, X1.*X3, X2.*X3, X1.*X1, X2.*X2, X3.*X3] ;

predictors = zscore(predictors) ;
text_vars = {'X1','X2','X3', 'X1*X2', 'X1*X3', 'X2*X3', 'X1*X1', 'X2*X2', 'X3*X3'} ;

results = cell(0,0) ;
idx=1 ;
daic=[] ;

for k = 1:size(predictors,2)
    vars = nchoosek(1:size(predictors,2),k) ;
    for i = 1:size(vars,1)
        if true
            var_val_M = [ones(numel(d_response),1),predictors(:,vars(i,:))] ;
            %precond = diag(1./sum(abs(var_val_M))) ;
            precond = eye(size(var_val_M,2)) ;
            lastwarn('') ;
            [b,bint,r,rint,stats] = regress(d_response, var_val_M*precond) ;
            [msgstr, msgid] = lastwarn ;
            if ~strcmp(msgid, '')
                warning_qty = warning_qty + 1 ;
            end
            b = precond*b ;
            real_k = k+1 ;
        else
            mdl=fitlm([XX(:,vars(i,:))], impact) ;
            tmp = anova(mdl, 'summary');
            tmp_F = tmp.F(isfinite(tmp.F));
            tmp_pValue = tmp.pValue(isfinite(tmp.pValue));
            stats = [mdl.Rsquared.Ordinary, tmp_F(1) , tmp_pValue(1), 0.0 ];
            r = sqrt(tmp.SumSq(3)) ;
            b = [] ;
            real_k = k+1 ;
        end
        n = numel(d_response) ;
        daic(idx) = 2*real_k + n*log(sum(r.*r)) +(2*real_k*(real_k+1))/(n-real_k-1);
        results(idx,:) = {stats, strjoin(text_vars(vars(i,:)), ', ') , b} ;
        idx=idx+1 ;
    end
end

true_daicc = daic - min(daic) ;
smallest_daicc = true_daicc <= Inf;

important_vars = {} ;
stats = [] ;
for i = find(smallest_daicc)
    R = results(i,:) ;
%    fprintf('Trophic level for graph: %d\n', R{1}) ;
%    fprintf('DAICc index=%d, value=%f\n', i, true_daicc(i)) ;
    stats(i,:) = R{1} ;
    
%    fprintf('R^2=%f, p=%f\n',stats(i,1), stats(i,3)) ;
%    fprintf('Coefficients: ')
    {true_daicc(i), R{2}}
%    fprintf('\n')
end

summary = table(true_daicc(smallest_daicc)',...
                stats(smallest_daicc,1),...
                stats(smallest_daicc,3),...
                results(smallest_daicc,2),...
                (1:length(smallest_daicc))') ;
summary.Properties.VariableNames = {'DAICc', 'R2', 'p', 'Imp_Vars', 'Index'} ;
sorted_summary = sortrows(summary,1);
sorted_summary(1:15,:)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table_all_pp = load_pp_files(experiment_pattern)
   file_list = dir(sprintf('%s/pp_*.mat',experiment_pattern)) ;
   file_table = struct2table(file_list) ;
   file_table = sortrows(file_table, {'folder', 'name'}) ;
   file_list = table2struct(file_table) ; % sort entries by folder+name
   
   table_all_pp.C = [] ;
   table_all_pp.C_measured = [] ;
   table_all_pp.S = [] ;
   table_all_pp.umig_rate = [] ;
   table_all_pp.sigma = [] ;
   table_all_pp.p_alpha = [] ;
   table_all_pp.p_beta_sub = [] ;
   table_all_pp.p_beta_div = [] ;
   table_all_pp.p_gamma = [] ;
   table_all_pp.b_alpha = [] ;
   table_all_pp.b_beta_sub = [] ;
   table_all_pp.b_beta_div = [] ;
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
       
       table_all_pp.C = [table_all_pp.C; all_C] ;
       table_all_pp.C_measured = [table_all_pp.C_measured; all_C_measured] ;
       table_all_pp.S = [table_all_pp.S; all_S] ;
       table_all_pp.umig_rate = [table_all_pp.umig_rate; all_umig_rate] ;
       table_all_pp.sigma = [table_all_pp.sigma; all_sigma] ;
       %   total_biomass = [table_all_pp.] ;
       table_all_pp.p_alpha = [table_all_pp.p_alpha; all_p_alpha] ;
       table_all_pp.p_beta_sub = [table_all_pp.p_beta_sub; all_p_beta_sub] ;
       table_all_pp.p_beta_div = [table_all_pp.p_beta_div; all_p_beta_div] ;
       table_all_pp.p_gamma = [table_all_pp.p_gamma; all_p_gamma] ;
       table_all_pp.b_alpha = [table_all_pp.b_alpha; all_b_alpha] ;
       table_all_pp.b_beta_sub = [table_all_pp.b_beta_sub; all_b_beta_sub] ;
       table_all_pp.b_beta_div = [table_all_pp.b_beta_div; all_b_beta_div] ;
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





















