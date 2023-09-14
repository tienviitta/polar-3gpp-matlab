function [b0o, b1o] = bintrl(b0i, b1i)
% BINTRL Bit Interleaver
%   b0o, b1o = BINTRL(b0i, b1i)
%   two vector input and two vector output 'perfect suffle' bit interleaver
%
%   'b0i' and 'b1i' are two input vectors of length N/2

bi = [b0i; b1i];
l_bi = length(bi);
bo = zeros(l_bi);

% Perfect Shuffle
for i = 0:l_bi-1
  if mod(i, 2)
    bo(i+1) = bi(floor(i/2)+floor(l_bi/2)+1);
  else
    bo(i+1) = bi(floor(i/2)+1);
  endif
endfor

% Vector Outputs
b0o = bo(1:floor(l_bi/2), 1);
b1o = bo(floor(l_bi/2)+1:l_bi, 1);

endfunction
