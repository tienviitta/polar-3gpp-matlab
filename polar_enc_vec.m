clc
clear all
close all

% Params
K = 8
N = 8

% Polar Encoding
%   3GPP TS 38.212 version 17.3.0 Release 17 ETSI TS 138 212 V17.3.0 (2022-09)
%   5.3.1.2 Polar encoding
bi = transpose(1:N);
bi(K+1:end) = 0
b0 = bi(1:N/2, 1);
b1 = bi(N/2+1:N, 1);
c0 = [1; 1; 1; 1]; % TODO: Compute control vectors!?

% Frozen Bit Insertion
for s = 1:log2(N)
  %   Butterfly Network
  [b0f, b1f] = bfly(b0, b1, c0);
  %   Bit Interleaver, e.g., 'perfect suffle'
  [b0, b1] = bintrl(b0f, b1f);
endfor
bo = [b0; b1]

