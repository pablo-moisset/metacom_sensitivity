function pm = possible_migs(patch, M, METACOM_SIM_P)
% matrix containing all possible migration delays for all possible
% migration events from/to patch
% pm is a sparse square matrix with patch_qty*richness rows (or columns)
% it is formatted and addresed just like mig_delay in metacom_sim()
% M is a structure containing the state of the metacommunity.
% Uses the same def as the matrix returned by create_metacom
% assumes patch is active

global MIG_R_TRS_IDX

%M.present 

if M.active(patch)~=1
assert(M.active(patch)==1) ;
end
%patch = double(patch) ;

patch_qty = length(M.landscape) ;
richness = length(M.A) ;

PR = patch_qty*richness ;
PR2 = patch_qty*richness^2 ;

pm = sparse(PR, PR) ;
%start from possible migrations into patch

for source_patch = 1:patch_qty
   if M.active(source_patch) && M.landscape(source_patch,patch)>0 %&& source_patch~=patch
       % compute species in source, but not in this patch
       % candidates is a logical vector
       candidates = (M.present(source_patch,:) - M.present(patch,:))==1 ;
       tmp = M.A ;
       tmp(~logical(M.present(patch,:)),:) = 0 ; %non present species are not food
       candidates = logical(((sum(tmp)>0) + M.basals).*candidates) ;

       assert(size(candidates,2)==richness) ;
       mr = M.unit_mig_rates./M.landscape(source_patch,patch) ;
%       candidates
%       mr
%source_patch
%patch
%       METACOM_SIM_P
%       METACOM_SIM_P.mig_rate_p(MIG_R_TRS_IDX)
       
       candidates = candidates & mr'>METACOM_SIM_P.mig_rate_p(MIG_R_TRS_IDX) ;
       if sum(candidates,2) > 0
          %These lines are ugly as heck
          %Transform source_patch, patch, candidates to linear indices in pm
          % idx is a vector because candidates is a vector
          idx = (patch-1)*PR2 + (source_patch-1)*richness + (find(candidates)-1)*(PR+1) + 1 ;
          pm(idx) = exprnd(M.landscape(source_patch,patch)./M.unit_mig_rates(candidates)', 1, sum(candidates,2)) ;
       end
   end
end

% Now we try to go from patch to anywhere else

for dest_patch = 1:patch_qty
   if M.active(dest_patch) && M.landscape(patch, dest_patch)
       % compute species in this patch but not in destination
       % candidates is a logical vector
       candidates = (M.present(patch,:) - M.present(dest_patch,:))==1 ;
       tmp = M.A ;
       tmp(~logical(M.present(dest_patch,:)),:) = 0 ; %non present species are not food
       candidates = logical(((sum(tmp)>0) + M.basals).*candidates) ;
       assert(size(candidates,2)==richness) ;

       mr = M.unit_mig_rates/M.landscape(patch,dest_patch) ;
       candidates = candidates & mr'>METACOM_SIM_P.mig_rate_p(MIG_R_TRS_IDX) ;

       if sum(candidates,2) > 0
%fprintf('from =%d, to=%d, candidates=', patch, dest_patch)
%fprintf('%d ', find(candidates))
%fprintf('\n');
          %These lines are ugly as heck
          %Transform patch, dest_patch, candidates to linear indices in pm
          % idx is a vector because candidates is a vector
          idx = (dest_patch-1)*PR2 + (patch-1)*richness + (find(candidates)-1)*(PR+1) + 1 ; 

          pm(idx) = pm(idx) + exprnd(M.landscape(patch,dest_patch)./M.unit_mig_rates(candidates)', 1, sum(candidates,2)) ;
       end
   end
end
%chg = full(pm)
end
