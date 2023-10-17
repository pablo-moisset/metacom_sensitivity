%UNK
%all_pa=[] ;
%all_ba=[] ;
%all_bg=[] ;
%all_mig_tsh = [] ;
%all_var = [] ;
%all_allo = [] ;


%mean_ba = [] ;
%se_ba   = [] ;
%mean_bg = [] ;
%se_bg   = [] ;

u_mig_tsh = unique(all_mig_tsh) ;
f_var = 0.4;

figure ;

legends = {} ;
for allo = unique(all_allo)'
   mean_pa = [] ;
   se_pa   = [] ;
   for mt = u_mig_tsh'
      idx = all_mig_tsh==mt & all_var==f_var & all_allo==allo;
      mean_pa = [mean_pa;mean(all_pa(idx))] ;
      se_pa = [se_pa; std(all_pa(idx))/sqrt(sum(idx))] ;
   end
   errorbar(u_mig_tsh, mean_pa, se_pa, '-o', 'MarkerSize',3.5, 'MarkerFaceColor', 'auto') ;
      hold on ;
  legends = {legends{:}, sprintf('$M_b=%g$', allo)} ; 
end
legend(legends,'interpreter','latex') ;
xlabel('Migration rate threshold', 'interpreter', 'latex') ;
ylabel('$P-\alpha$', 'interpreter', 'latex') ;

