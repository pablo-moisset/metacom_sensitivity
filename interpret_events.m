function [M,colonizations,extinctions] = interpret_events(E, metacom)
% E is a cell array of events. metacom is the initial state of the metacommunity
% This function returns a cell array of metacommunity states after each event
   
   % Preallocate
   M = cell(length(E)+1,1) ;
   colonizations = zeros(length(E)+1,1) ;
   extinctions   = zeros(length(E)+1,1) ;
   
   M{1} = metacom ;
   
   for i = 1:length(E)
       e = E{i} ;
       metacom.timestamp = e.timestamp ;
       metacom.last_touched_patch = e.last_touched_patch ;
       metacom.last_event_type = e.last_event_type ;
       if e.last_event_type == 1 %Activate event
           metacom.active(e.last_touched_patch) = true ;
           %metacom.present(e.last_touched_patch,:) = false ;
       elseif e.last_event_type == 2  %dectivate event
           metacom.active(e.last_touched_patch) = false ;
           metacom.present(e.last_touched_patch,:) = false ;
       elseif e.last_event_type == 3  %migration attemp event
           metacom.present(e.last_touched_patch,e.added) = true ;
           metacom.present(e.last_touched_patch,e.removed) = false ;
           colonizations(i) = sum(e.added) ;
           extinctions(i) = sum(e.removed) ;
       else
           assert(false) ;
       end
       M{i+1} = metacom ; %Attach new metacommunity state
   end
end
