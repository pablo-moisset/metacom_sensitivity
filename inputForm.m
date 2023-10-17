function s = inputForm(A, f)
%s = inputForm(A, name)
% Create assignment code from disp output
%
% Inputs
%   A: anything disp can handle
%   name: variable A is assigned to
% Output
%   s:  string with expression 'name = A;'

fprintf(f,"[") ; 
for i=1:size(A,1)
   fprintf(f,"%d ", A(i,:)) ; 
   fprintf(f,";...\n") ; 
end     
fprintf(f,"]\n") ; 
