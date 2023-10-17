function metaplot(topology, response_vars, filter_variability)
%topology='le20v_d-0.2T1.5'
%topology='le20v_d-0.2T0'
%topology='le20v_d-2T1.5'
%topology='le20v_d-2T0'
%lc_pattern='0.3'

%response_vars = {'p_alpha'} ;
%response_vars = {'p_beta_sub'} ;
%response_vars = {'p_beta_div'} ;
%response_vars = {'p_gamma'} ;
%response_vars = {'shannon_alpha'}
%response_vars = 'shannon_beta'
%response_vars = {'shannon_gamma'}
%response_vars = 'b_alpha'
%response_vars = 'b_beta_sub'
%response_vars = 'b_beta_div'
%response_vars = {'b_gamma'} ;

%ivar = 'Wet' ;
%ivar = 'Var' ;
ivar = 'Mig' ;

%panel_n = 'Wet' ;
panel_n = 'Var' ;
%panel_n = 'Mig' ;

%filter_function = @(wetness, variability, mig_rate, ext_tsh) variability==0.0 & ext_tsh==0.2;
%filter_function = @(wetness, variability, mig_rate, ext_tsh) true(size(variability));
filter_function = @(wetness, variability, mig_rate, ext_tsh) abs(variability-filter_variability)<1e-5;
%filter_function = @(wetness, variability, mig_rate, x_rate) mig_rate==30 ;

all_pp = load_pp_files(topology, filter_function) ;

for i=1:numel(response_vars)
    response = response_vars{i}
 %   close all ;
    f_handles=[] ;
    p_handles=[] ;

    panel_n_values = find_panel_n_values(all_pp, panel_n) ;
    
    for j = 1:length(panel_n_values)
       [response_data, ivar_data, complexity] = filter_panel_values(all_pp, response, ivar, panel_n, panel_n_values(j)) ; 
       [ivar_data,I] = sort(ivar_data) ;
       response_data = response_data(I) ;
       response_data = reshape(response_data,[],numel(unique(ivar_data)));
       errorbar(log10(unique(ivar_data)),mean(response_data),std(response_data)/sqrt(size(response_data,2))) ;
       %plotbanda_mean(log10(unique(ivar_data))',mean(response_data),std(response_data)/sqrt(size(response_data,2)),'xLabel','yLabel','Title')
       title(topology, 'Interpreter', 'none')
       return
       assert(false) ;
       [fh ph]=make_plot(topology, response_data, ivar_data, complexity, response, ivar, panel_n, panel_n_values(j)) ;
       f_handles = [f_handles, fh] ;
       p_handles = [p_handles, ph] ;
    end

   % Create combined figure

   fig_name = sprintf('T:%s, response: %s', topology, response) ;
   hf_main = figure('Name', fig_name, 'NumberTitle','off');
   npanels = numel(f_handles);
   hp_sub = nan(1,npanels);
   % Copy over the panels
   for idx = 1:npanels
      hp_sub(idx) = copyobj(p_handles(idx),hf_main);
      set(hp_sub(idx),'Position',[0,(npanels-idx)/npanels,1,1/npanels]);
   end
   fig_name = sprintf('T%s_resp_%s.fig', topology, response) ;
   savefig(hf_main, fig_name) ;   
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table_all_pp = load_pp_files(topology, filter_function)
   file_list = dir(sprintf('outputs_%s/pp_*.mat',topology)) ;

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
   end
   rows_to_keep = filter_function(table_all_pp.wetness,...
       table_all_pp.variability,...
       table_all_pp.umig_rate,... 
       table_all_pp.ext_tsh) ;
   table_all_pp.C = table_all_pp.C(rows_to_keep) ;
   table_all_pp.C_measured = table_all_pp.C_measured(rows_to_keep);
   table_all_pp.S = table_all_pp.S(rows_to_keep) ;
   table_all_pp.umig_rate = table_all_pp.umig_rate(rows_to_keep) ;
   table_all_pp.sigma = table_all_pp.sigma(rows_to_keep) ; 
   %   total_biomass = [table_all_pp.] ;
   table_all_pp.p_alpha = table_all_pp.p_alpha(rows_to_keep) ;
   table_all_pp.p_beta_sub = table_all_pp.p_beta_sub(rows_to_keep) ;
   table_all_pp.p_beta_div = table_all_pp.p_beta_div(rows_to_keep) ;
   table_all_pp.p_gamma = table_all_pp.p_gamma(rows_to_keep) ;
   table_all_pp.b_alpha = table_all_pp.b_alpha(rows_to_keep) ;
   table_all_pp.b_beta_sub = table_all_pp.b_beta_sub(rows_to_keep) ; 
   table_all_pp.b_beta_div = table_all_pp.b_beta_div(rows_to_keep) ; 
   table_all_pp.b_gamma = table_all_pp.b_gamma(rows_to_keep) ;
   table_all_pp.ext_tsh = table_all_pp.ext_tsh(rows_to_keep) ;
   table_all_pp.variability = table_all_pp.variability(rows_to_keep) ;
   table_all_pp.shannon_alpha = table_all_pp.shannon_alpha(rows_to_keep); 
   table_all_pp.shannon_beta = table_all_pp.shannon_beta(rows_to_keep) ;
   table_all_pp.shannon_gamma = table_all_pp.shannon_gamma(rows_to_keep) ;
   table_all_pp.wetness = table_all_pp.wetness(rows_to_keep) ;
end

function panel_n_values = find_panel_n_values(all_pp, panel_n)
   if strcmp(panel_n, 'Var')
       tmp = all_pp.variability ;
   elseif strcmp(panel_n, 'Mig')
       tmp = all_pp.mig_rate ;
   elseif strcmp(panel_n, 'Wet')
       tmp = cellfun(@(x)x(2), all_pp.wetness) ;
   else
       assert(false) ;
   end
   panel_n_values = unique(tmp) ;
end

function [response_data, ivar_data, complexity] = filter_panel_values(all_pp, response, ivar, panel_n, value)
   if strcmp(panel_n, 'Var')
       rows_to_keep = all_pp.variability == value ;
   elseif strcmp(panel_n, 'Mig')
       rows_to_keep = all_pp.mig_rate == value ;
   elseif strcmp(panel_n, 'Wet')
       rows_to_keep = cellfun(@(x)x(2), all_pp.wetness) == value  ;
   else
       assert(false) ;
   end

   %response_vars = {'p_alpha','p_beta','p_gamma','shannon_alpha','shannon_beta','shannon_gamma','b_alpha','b_cv', 'b_gamma'} 

   command_str = sprintf('response_data = all_pp.%s(rows_to_keep);', response) ;
   eval(command_str) ;
   
   if strcmp(ivar, 'Var')
       ivar_data = all_pp.variability(rows_to_keep) ;
   elseif strcmp(ivar, 'Mig')
       ivar_data = all_pp.umig_rate(rows_to_keep) ;
   elseif strcmp(ivar, 'Wet')
       tmp = all_pp.wetness(rows_to_keep) ;
       ivar_data = 1-cellfun(@(x)x(2), tmp) ;
   else
       assert(false) ;
   end
   
   fC = all_pp.C_measured(rows_to_keep) ;
   fS = all_pp.S(rows_to_keep) ;
   fsigma = all_pp.sigma(rows_to_keep) ;
   complexity = fsigma.*sqrt(fC.*fS) ;
end
