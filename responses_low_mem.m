function [patch_present, persistence_gamma, t] = responses_low_mem(T, include_sources)
% T is a 1d cell array of metacom structures
% present_qty(p,k):  number of species present in site p at time index k
% persitence_gamma(k):  number of species present in the metacommunity at time index k
% This implementation is ugly and computes the response variables block
% by block to save memory

richness_master = length(T{1}.A) ; 

T  = T(1:end-1) ; %remove the last element to avoid inf timestamp
CM = T{1}.CM ;
R  = T{1}.R ;

L = length(T) ;
patch_qty = size(T{1}.present,1) ;

%consider_patches = 1:patch_qty ;
%if ~include_sources
%   consider_patches = setdiff(consider_patches, T{1,1}.sources) ;
%end

fprintf('Nbr of events = %d\n', L) ;

%tspan = T{end}.timestamp - T{1}.timestamp ;
t = zeros(L,1) ;
patch_present = zeros(patch_qty,L) ;
persistence_gamma = zeros(L,1) ;

for k = 1:L %k is a time index
    t(k) = T{k}.timestamp ;
    P = T{k}.present ; %P(i,j)==1 iff species j is in patch i at time t(k) 
    patch_present(:,k) = sum(P,2)' ; % Richnes at time t(k), for every patch
    tmp=P(2:end,:) ;
    sum(tmp(:));
    persistence_gamma(k) = sum(sum(P(2:end,:))>0) ;
end
%persistence_gamma
end
