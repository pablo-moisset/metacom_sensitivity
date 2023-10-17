topologies={...
'p50le100S30c5x1u_ver3_T0.001_D-0.2', 'p50le100S30c5x1u_ver3_T0.35_D-0.2', 'p50le100S30c5x1u_ver3_T0.7_D-0.2',...
'p50le100S30c5x1u_ver3_T0.001_D-0.6', 'p50le100S30c5x1u_ver3_T0.35_D-0.6', 'p50le100S30c5x1u_ver3_T0.7_D-0.6',...
'p50le100S30c5x1u_ver3_T0.001_D-1', 'p50le100S30c5x1u_ver3_T0.35_D-1', 'p50le100S30c5x1u_ver3_T0.7_D-1'...
}

%lc_pattern='0.3'

%response_vars = {'p_alpha','p_beta_sub','p_beta_div','p_gamma','shannon_alpha','shannon_beta','shannon_gamma','b_alpha','b_beta_sub','b_beta_div', 'b_gamma'} 
r_var = {'p_alpha'} ;

variabilities = linspace(0,0.8,5) ;
%variabilities = linspace(0,0.5,10) ;

figure;
%h = figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
%ax = subplot(1,1,1);
%ax.Position = [0.25 0.25 0.65 0.65];
%ax.ActivePositionProperty = 'position';
%ax.TickLabelInterpreter='latex';

for tpl_idx = 1:numel(topologies)
   tpl = topologies(tpl_idx);
   subplot(3,3,tpl_idx) ;
   metaplot(tpl{1}, r_var, variabilities(1)) ;
   xlabel('log10(mig rate a)')
   ylabel(r_var, 'Interpreter', 'none')
   hold on ;
   for v = variabilities
       metaplot(tpl{1}, r_var, v) ;
   end
end