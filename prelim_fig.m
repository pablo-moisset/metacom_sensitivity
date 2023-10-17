

% Regular expression identifies tseries files with the same migration rate
%table0 = load_pp_files('MOD_landscape/S45_p4/outputs_p50le100S45c5x50u_var0_T0.6_D-0.6_px2_CON0.2',   '.*18085.[5-9]_.*');
%table1 = load_pp_files('MOD_landscape/S45_p4/outputs_p50le100S45c5x50u_var0.5_T0.6_D-0.6_px2_CON0.2', '.*17938.[5-9]_.*');
table0 = load_pp_files('MOD_landscape/S45_p4/outputs_p100le200S45c5x50u_var0_T0.2_D-0.467_px2_CON0.2',   '.*15394.[0-4]_.*');
table1 = load_pp_files('MOD_landscape/S45_p4/outputs_p100le200S45c5x50u_var0.5_T0.2_D-0.467_px2_CON0.2', '.*15247.[0-4]_.*');
%table0 = load_pp_files('MOD_landscape/S45_p4/outputs_p10le20S45c5x50u_var0_T0.2_D-0.467_px2_CON0.2',   '.*16276.[0-4]_.*');
%table1 = load_pp_files('MOD_landscape/S45_p4/outputs_p10le20S45c5x50u_var0.5_T0.2_D-0.467_px2_CON0.2', '.*16129.[0-4]_.*');
%5-9
%15-19
rvar0 = extract_rv(table0) ; %Convert from table fromat (as returned by load_pp_files
rvar1 = extract_rv(table1) ; %To a more convenient matrix format
% Each column of rvar0 and rvar1 are associated with a response var
% while the rows have the values

Z = rvar1./rvar0 ;
mean_vals = [100*(mean(Z)-1)] ;
se_vals = [100*sqrt((mean(rvar1.^2).*mean(1./rvar0.^2)-mean(rvar1).^2.*mean(1./rvar0).^2) / size(rvar1,1) ) ];  

%table0 = load_pp_files('MOD_landscape/S45_p4/outputs_p50le100S45c5x50u_var0_T0.6_D-0.6_px2_CON0.2',   '.*18085.1[5-9]_.*');
%table1 = load_pp_files('MOD_landscape/S45_p4/outputs_p50le100S45c5x50u_var0.5_T0.6_D-0.6_px2_CON0.2', '.*17938.1[5-9]_.*');
table0 = load_pp_files('MOD_landscape/S45_p4/outputs_p100le200S45c5x50u_var0_T0.2_D-0.467_px2_CON0.2',   '.*15394.2[0-4]_.*');
table1 = load_pp_files('MOD_landscape/S45_p4/outputs_p100le200S45c5x50u_var0.5_T0.2_D-0.467_px2_CON0.2', '.*15247.2[0-4]_.*');
%table0 = load_pp_files('MOD_landscape/S45_p4/outputs_p10le20S45c5x50u_var0_T0.2_D-0.467_px2_CON0.2',   '.*16276.2[0-4]_.*');
%table1 = load_pp_files('MOD_landscape/S45_p4/outputs_p10le20S45c5x50u_var0.5_T0.2_D-0.467_px2_CON0.2', '.*16129.2[0-4]_.*');
%5-9
%15-19
rvar0 = extract_rv(table0) ;
rvar1 = extract_rv(table1) ;

Z = rvar1./rvar0 ;
mean_vals = [mean_vals; 100*(mean(Z)-1)] ;
se_vals = [se_vals;...
100*sqrt((mean(rvar1.^2).*mean(1./rvar0.^2)-mean(rvar1).^2.*mean(1./rvar0).^2) / size(rvar1,1) ) ];  

figure
bh = bar(mean_vals', 'BarWidth', 1) ;
set(bh(1),'facecolor',[100 170 166]/255.0);
set(bh(2),'facecolor',[255 194 0]/255.0);
bh(1).EdgeColor = 'none' ;
bh(2).EdgeColor = 'none' ;
hold on
ngroups = size(mean_vals', 1);
nbars = size(mean_vals', 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    h = errorbar(x, mean_vals(i,:), se_vals(i,:));
    h.LineStyle = 'none' ;
    h.Color = 'Black' ;
end
hAxes = get(gca, 'XAxis');
hAxes.TickLabelInterpreter = 'latex';
hAxes.FontSize = 13;
xticklabels({'$\mathcal{P}_\alpha$','$\mathcal{P}_\beta$', '$\mathcal{P}_\gamma$',...
             '$\mathcal{B}_\beta$', '$\mathcal{B}_\gamma$'})
%             '$\mathcal{H}_\alpha$','$\mathcal{H}_\beta$', '$\mathcal{H}_\gamma$',...
%ylim([-200 1000]) ;
xlim([0,6]) ;
xlabel('Response variable', 'fontsize',12)
ylabel('Relative change (%)', 'fontsize',12)

legend({'$a=30$', '$a=3000$'}, 'interpreter', 'latex')
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract response variables, i.e. go from table to matrix format
function r = extract_rv(table)
r=[table.p_alpha, table.p_beta, table.p_gamma, ...
   table.b_beta, table.b_gamma] ;
%   table.shannon_alpha, table.shannon_beta, table.shannon_gamma,...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table_all_pp = load_pp_files(experiment_pattern,expression)
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
    if ~isempty(regexp(file_list(i).name, expression))
        fprintf('File nbr. %d: %s\n',i,file_list(i).name)
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
end
