function status = tvwrite(fn, tv)
% TVWRITE Write vector to a testvector file
%   status = tvwrite(fn, tv)

status = 1;

% Open file
fid = fopen(fn, "w");

% Write scalar vector to a file line by line
l_tv = length(tv);
for i = 1:l_tv
  fprintf(fid, "%d\n", tv(i));
endfor

% Close file
fclose(fid);

status = 0;

endfunction

