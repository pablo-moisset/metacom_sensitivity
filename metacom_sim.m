% Fixes a bug in prev sims
function E = metacom_sim(M, Q, METACOM_SIM_P)
% Runs the event driven simulation until running out of events.
% M is a structure representing the metacommunity (initial state and parameters)
% Q is a list of patch on, off events

   MIGRATION_EVENT_TYPE = 3 ;

   M.active = M.active*0 ; %hack, source patches are activated as events in Q
                           % if they are "preactivated" it causes problems  
                           %an alternative would have been to remove source
                           %activation events from q
   E = {} ;
   E_len = 0 ;

   richness = length(M.A) ;
   patch_qty = length(M.landscape) ;

   mig_delay = sparse(richness*patch_qty, richness*patch_qty) ;
   t = 0 ;
   last_t = last_NI_patch_event_time(Q) ;

   while t < last_t % <= maybe ?????
      nz_idx = find(mig_delay) ;
      [delay, idx] = min(mig_delay(nz_idx)) ;
      mig_idx = nz_idx(idx) ;
      if isempty(delay) %empty queue of migrations
          mig_t = Inf ;
      else
          mig_t = t + delay ; % Absolute time for earliest migration event
      end

      first_on_off = Q{1} ;

      if first_on_off.time <= mig_t % patch on/off event is first
          % Correct mig_delays, because we are advancing the time by (first_on_off.time - t)
          mig_delay(nz_idx) = mig_delay(nz_idx) - (first_on_off.time - t) ;
          p = first_on_off.patch ; % Index for patch going on or off

         if first_on_off.type == 1
             M.active(p) = true; 
             %COMPUTE POSSIBLE MIGRATIONS into this patch
             pm = possible_migs(p, M, METACOM_SIM_P) ;
%             full(mig_delay)
             assert( sum(sum(mig_delay.*pm))==0 ) ;
             mig_delay = mig_delay + pm ;
         elseif first_on_off.type == 2
            p = first_on_off.patch ;  % Index for patch going off
            M.active(p) = false ;
            M.present(p, :) = false ;
            mig_delay( (p-1)*richness + (1:richness) , : ) = 0 ; % Kills all events
            mig_delay( : , (p-1)*richness + (1:richness) ) = 0 ; % of migrations to/from patch p 
         else
            assert(false) ;
         end
         E_len = E_len+1 ;

         %new_sim_event = sim_event ;
         new_sim_event = struct() ;
         new_sim_event.timestamp = first_on_off.time ;
         new_sim_event.last_event_type = first_on_off.type ;
         new_sim_event.last_touched_patch = p ;

         E{E_len} = new_sim_event ; % E{E_len} ;
         t = first_on_off.time ;
         Q = Q(2:end) ;
      else % migration event
         % Correct mig_delays, because we are advancing the time by (first_on_off.time - t)
         mig_delay(nz_idx) = mig_delay(nz_idx) - delay ;
         mig_delay(mig_idx) = 0 ; % probably redundant because of the previous line
         [row,col] = ind2sub(size(mig_delay), mig_idx) ;
         from_idx = floor((row-1)/richness)+1 ;
         to_idx   = floor((col-1)/richness)+1 ;
         assert(M.active(from_idx)==1 & M.active(to_idx)==1) ;
         species_idx = mod(uint64(row-1), uint64(richness)) + 1 ;
         assert ( species_idx == mod(uint64(col-1), richness)+1 ) ;

         [to_add, to_remove] = eq_and_try_add(to_idx, species_idx, M,...
                                              METACOM_SIM_P.check_feasibility) ;
         assert(size(to_add,2)==richness) ;
         assert(size(to_remove,2)==richness) ;
         M.present(to_idx, to_add) = true ;
         M.present(to_idx, to_remove) = false ;

         if sum(to_add,2)>0 || sum(to_remove,2)>0
            pm = possible_migs(to_idx, M, METACOM_SIM_P) ;
            mig_delay(double(to_idx-1)*richness+(1:richness),:) = 0 ;
            mig_delay(:,double(to_idx-1)*richness+(1:richness)) = 0 ;
            
%            assert( all(mig_delay(:).*pm(:)==0) ) ;
            mig_delay = mig_delay + pm ;

            E_len = E_len+1 ;

            %new_mig_event = mig_event ;
            new_mig_event = struct() ;
            new_mig_event.timestamp = mig_t ; 
            new_mig_event.last_event_type = MIGRATION_EVENT_TYPE ;
            new_mig_event.last_touched_patch = to_idx ;
            new_mig_event.added = to_add ;
            new_mig_event.removed = to_remove ;

            E{E_len} = new_mig_event ;
         end
         t = mig_t ;
      end
   end
%   E
end

function t = last_NI_patch_event_time(Q)
% Largest finite timestamp among all on-off events
   t = 0 ;
   for i=1:numel(Q)
       t_event = Q{i}.time ;
       if t_event~=Inf && t_event > t
           t = t_event ;
       end
   end
end

