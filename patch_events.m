function [events, event_matrix] = patch_events(n, permanent, T, years, length_on, variability)
% List of on/off events for patches represented as a cell array
% the list is sorted by time
% each entry is a structure containing
%   type=1 or type=2 representing on/off event respectively
%   time: timestamp of event
%   patch: the patch index number, an integer between 1 and n
% Parameters
%    n: number of patches (positive integer)
%    permanent: 1d array containing patch index numbers indicating those
%               patches are always on.
%    T: length of one year
%    years: number of years to generate
%    length_on: length of the "on" period, without considering overlaps
%    variability: degree of asynchrony among patches

assert( n > 0 ) ;
assert( all(permanent>=1) && all(permanent<=n) ) ;
assert( T > 0 ) ;
assert( years >= 1 ) ;
assert( variability >= 0 ) ;

event_matrix = ones(1, 2*years) ;
for i = 1:years
    event_matrix(2*i-1) = i ;
    event_matrix(2*i) = i ;
end

event_matrix = repmat(event_matrix, n, 1)*T;

for i = 1:years
    displacement = (2*rand(n,1)-1.0)*variability ;
    event_matrix(:,2*i-1) = event_matrix(:,2*i-1) + displacement ;
    event_matrix(:,2*i) = event_matrix(:,2*i) + length_on + displacement ; 

    event_matrix(permanent, 2*i-1) = Inf ;
    event_matrix(permanent, 2*i) = Inf ;
end

event_matrix = event_matrix - min(min(event_matrix)) ;

for i = 1:years
    event_matrix(permanent, 2*i-1) = 0 ;
end

events = {} ;
j=0 ;
event_timestamp = [] ;

for k = 1:n
    [timestamp, I] = sort(event_matrix(k,:), 2) ;
    type = repmat([1 -1], 1, years) ;
    tmp = cumsum(type(I), 2) ;

    for i = 1:2*years
        if tmp(i)==1 && (i==1 || tmp(i-1)==0)
            start_on = i ;
        elseif tmp(i)==0
            e.patch = k ;
            e.type = 1 ;
            e.time = timestamp(start_on) ;
            j = j + 1 ; events{j}= e ;
            event_timestamp = [event_timestamp, e.time] ;
            e.type = 2 ;
            e.time = timestamp(i) ;
            j = j + 1 ; events{j}= e ;
            event_timestamp = [event_timestamp, e.time] ;
%            sprintf("Adding from %f to %f", start_on, i)
        end
    end
end
[s,I] = sort(event_timestamp, 2);
%s
%events
events = events(I);
end