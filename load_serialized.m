function result = load_serialized(filename_in)
   assert(exist(filename_in, 'file') == 2) ;
   tmpName = tempname('/tmp') ;
   r = system(sprintf('gunzip -c %s >%s', filename_in, tmpName)) ;
   [f,errmsg] = fopen(tmpName, 'r') ;
   bin_data = fread(f, 'uint8') ;
   result = hlp_deserialize(bin_data) ;
   clear bin_data ;
   fclose(f) ;
   delete(tmpName) ;
%   assert(r==0 && errmsg=="" ) ;
end
