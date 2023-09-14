function [b0o, b1o] = bfly(b0i, b1i, cs)
% BFLY Butterfly Network
%   b0o, b1o = BFLY(b0i, b1i, cs)
%   two vector input and two vector output butterfly network
%
%   'b0i' and 'b1i' are two input vectors of length N/2 and 'cs' is the control
%   for the butterfly network of length N/2

b0o = zeros(length(cs), 1);
b1o = zeros(length(cs), 1);

for i = 1:length(cs)
  if cs(i) == 1
    b0o(i) = b0i(i);
  else
    b0o(i) = b1i(i);
  endif
  if cs(i) == 1
    b1o(i) = b1i(i);
  else
    b1o(i) = b0i(i);
  endif
endfor

endfunction
