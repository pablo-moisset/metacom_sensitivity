function M = create_metacom(landscape, sources, non_intermitent, A, mig_rate_p)
% Metacommunity describing species in each patch at the initial time
% M.landscape: distance matrix among patches.
%    M.landscape(i,j)~=0 iff from i we can reach j in one hop.
% M.A: adjacency 0-1 master matrix. M.A(i,j)==1 iff j preys on i
% M.sources: vector of integers representing (as indices) the subset of 
%      patches that act as sources (A.K.A. continents, where all species are always present)
% M.non_intermitent: vector of integers representing (as indices) the subset of 
%      patches that never dry out (are not continents)
% M.unit_mig_rates: vector representing migration rates per each species (for distance 1).
% M.present: Matrix representing present species in each patch.
%    M.present(p,i) == true iff i is present in patch p, and is false otherwise;
% M.timestamp: set to default time zero
% M.active: boolean vector describing the set of patches currently at state "on"
%    We set all entries corresponding to sources to 1.
% M.basals: 1d 0-1 vectors flagging basal species

assert(all(sources>=1) && all(sources<=size(landscape,1))) ;
assert(all(non_intermitent>=1) && all(non_intermitent<=size(landscape,1))) ;
assert(isempty(intersect(sources,non_intermitent))) ;

assert(length(mig_rate_p)==3) ;

global MIG_R_KM0_IDX MIG_R_BA_IDX MIG_R_TRS_IDX

M.A = A ;
M.landscape = landscape ;
M.sources = sources ;
M.non_intermitent = non_intermitent ;
M.present = false(length(landscape), size(A,1)) ;
M.present(sources,:) = true ;
M.timestamp = 0 ;
M.active = zeros(size(landscape,1),1) ;
M.active(sources) = 1 ;
%M.active(non_intermitent) = 1 ;

M.basals = (sum(A)-diag(A)')==0 ;

M.last_event_type = NaN ;
M.last_touched_patch = NaN ;

% Compute migration rates using allometry
TP = TP_shortestpath(A) ;
M.unit_mig_rates = mig_rate_p(MIG_R_KM0_IDX)*mig_rate_p(MIG_R_BA_IDX).^TP ;
max_mr = max(M.unit_mig_rates) ;

nz = M.landscape > 0 ;
tmp = M.landscape(nz) ;
tmp(max_mr./tmp<mig_rate_p(MIG_R_TRS_IDX)) = 0 ;

M.landscape(nz) = tmp;
end
