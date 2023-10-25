function f = PDCCH_encoder_ref(a, E, RNTI, tc)
% PDCCH_ENCODER Polar encoder for the Physical Downlink Control Channel (PDCCH) of 3GPP New
% Radio, as defined in Section 7.3 of TS38.212. Implements the zero-
% padding to increase the length of short payloads to 12 bits of Section 7.3.1,
% the Cyclic Redudancy Check (CRC) attachment of Section 7.3.2, the channel
% coding of Section 7.3.3 and the rate matching of Section 7.3.4. Note that
% this code does not implement the DCI bit sequence generation of Section
% 7.3.1.
%   f = PDCCH_ENCODER(a, E) encodes the information bit sequence a, in
%   order to obtain the encoded bit sequence e.
%
%   a should be a binary row vector comprising A number of bits, each
%   having the value 0 or 1. A should be in the range 1 to 140. The first
%   input bit corresponds to a_0 from Section 7.3.1 of TS38.212, while the
%   last input bit corresponds to a_A-1.
%
%   E should be an integer scalar. It specifies the number of bits in the
%   encoded bit sequence, where E should greater than A.
%
%   RNTI should be a binary row vector comprising 16 bits, each having the
%   value 0 or 1. If this parameter is omitted, then ones(1,16) will be
%   used for the RNTI. The first bit corresponds to x_rnti,0 from Section
%   7.3.2 of TS38.212, while the last bit corresponds to x_rnti,15.
%
%   f will be a binary row vector comprising E number of bits, each having
%   the value 0 or 1. The first output bit corresponds to f_0 from Section
%   7.3.4 of TS38.212, while the last output bit corresponds to
%   f_E-1.
%
%   See also PDCCH_DECODER
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you
% can redistribute it and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version. This
% program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.

addpath 'components'

if nargin == 2
    RNTI = ones(1,16);
end
if length(RNTI) ~= 16
    error('RNTI length should be 16');
end



A = length(a);

if A == 0
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no less than 1.');
end
if A > 140
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no greater than 140.');
end
if E > 8192
    error('polar_3gpp_matlab:UnsupportedBlockLength','E should be no greater than 8192.');
end


% The CRC polynomial used in 3GPP PBCH and PDCCH channel is
% D^24 + D^23 + D^21 + D^20 + D^17 + D^15 + D^13 + D^12 + D^8 + D^4 + D^2 + D + 1
crc_polynomial_pattern = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
P = length(crc_polynomial_pattern)-1;

% If a contains fewer than 12 bits, increase its length to 12 by padding it
% with zeros
if A < 12
    a = [a,zeros(1,12-length(a))];
    K = 12+P;
else
    K = A+P;
end

% Determine the number of bits used at the input and output of the polar
% encoder kernal.
N = get_3GPP_N(K,E,9); % n_max = 9 is used in PBCH and PDCCH channels

% Get the 3GPP CRC interleaver pattern.
crc_interleaver_pattern = get_3GPP_crc_interleaver_pattern(K);

% Get the 3GPP rate matching pattern.
[rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K,N,E);

% Get the 3GPP sequence pattern.
Q_N = get_3GPP_sequence_pattern(N);

% Get the 3GPP information bit pattern.
info_bit_pattern = get_3GPP_info_bit_pattern(K, Q_N, rate_matching_pattern, mode);

% Perform Distributed-CRC-Aided polar encoding.
%f = DS1CA_polar_encoder(a,crc_polynomial_pattern, RNTI, crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern);
%A = length(a);
%P = length(crc_polynomial_pattern)-1;
%K = length(crc_interleaver_pattern);
%N = length(info_bit_pattern);

%if A+P ~= K
%    error('A+P should equal K');
%end
%if log2(N) ~= round(log2(N))
%    error('N should be a power of 2');
%end
%if sum(info_bit_pattern) ~= K
%    error('info_bit_pattern should contain K number of ones.');
%end
%if max(rate_matching_pattern) > N
%    error('rate_matching_pattern is not compatible with N');
%end
%if P < length(crc_scrambling_pattern)
%    error('polar_3gpp_matlab:UnsupportedBlockLength','K should be no less than the length of the scrambing pattern');
%end

% Generate the CRC bits.
G_P = get_crc_generator_matrix(A+P,crc_polynomial_pattern);
crc_bits = mod([ones(1,P),a]*G_P,2);

% Scramble the CRC bits.
scrambled_crc_bits = xor(crc_bits,[zeros(1,P-length(RNTI)),RNTI]);

% Append the scrambled CRC bits to the information bits.
b = [a, scrambled_crc_bits];

% Interleave the information and CRC bits.
c = b(crc_interleaver_pattern);

% Position the information and CRC bits within the input to the polar
% encoder kernal.
u = zeros(1,N);
u(info_bit_pattern) = c;

% Perform the polar encoder kernal operation.
G_N = get_G_N(N);
d = mod(u*G_N,2);

% Extract the encoded bits from the output of the polar encoder kernal.
f = d(rate_matching_pattern);

% Testvectors
printf("Distributed-CRC-Aided Polar encoding: A:%d, P:%d, K:%d, E:%d, N:%d, mode:%s\n", A, P, K, E, N, mode);
tvwrite([tc "/" "params.txt"], [A, P, K, E, N]);
tvwrite([tc "/" "info_bits.txt"], a);
tvwrite([tc "/" "crc_gen_m.txt"], G_P(:));
tvwrite([tc "/" "enc_gen_m.txt"], G_N(:));
tvwrite([tc "/" "rnti_bits.txt"], RNTI);
tvwrite([tc "/" "crc_interleaver_pattern.txt"], crc_interleaver_pattern-1);
tvwrite([tc "/" "rate_matching_pattern.txt"], rate_matching_pattern-1);
tvwrite([tc "/" "info_bit_pattern.txt"], info_bit_pattern);
tvwrite([tc "/" "crc_bits.txt"], crc_bits);
tvwrite([tc "/" "scrambled_crc_bits.txt"], scrambled_crc_bits);
tvwrite([tc "/" "info_crc_bits.txt"], b);
tvwrite([tc "/" "info_intrl_bits.txt"], c);
tvwrite([tc "/" "frozen_bits.txt"], u);
tvwrite([tc "/" "enc_bits.txt"], d);
tvwrite([tc "/" "rm_bits.txt"], f);

