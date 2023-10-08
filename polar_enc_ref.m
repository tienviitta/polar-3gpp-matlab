function status = polar_enc_ref(A, E, tc, seed)
% Polar encoder for the Physical Uplink Control Channel (PUCCH) and the
% Physical Uplink Shared Channel (PUSCH) of 3GPP New Radio, as defined in
% Section 6.3 of TS38.212. Implements the code block segmentation and
% Cyclic Redudancy Check (CRC) attachment of Sections 6.3.1.2.1 and 6.3.2.2.1,
% the channel coding of Sections 6.3.1.3.1 and 6.3.2.3.2, the rate matching of
% Sections 6.3.1.4.1 and 6.3.2.4.1, as well as the code block concatenation of
% Sections 6.3.1.5.1 and 6.3.2.5.1. Note that this code does not implement the
% UCI bit sequence generation of Sections 6.3.1.1 and 6.3.2.1, the
% determination of the encoded block length E_UCI of Sections 6.3.1.4.1 and
% 6.3.2.4.1, or the multiplexing of Sections 6.3.1.6 and 6.3.2.6. Also, this
% code does not implement the small block lengths, which are detailed in
% Sections 6.3.1.2.2, 6.3.1.3.2, 6.3.1.4.2, 6.3.2.2.2, 6.3.2.3.2 and 6.3.2.4.2.
%   POLAR_ENC_REF(A, E, TC, SEED) generates testvectors for Polar encoding DSP SW
%   implementation.
%
%   A should be in the range 12 to 1706. E should be an integer scalar. It
%   specifies the number of bits in the encoded bit sequence, where E should be
%   no greater than 8192 if A<360 and no greater than 16384 if A>=360.
%
%   Testcase tc is a string to use as a directory for the testvectors.
%

status = 1;

% Seed the random number generator
rng(seed);

% Generate a random block of bits
a = round(rand(1,A));

%%% Polar Encoding
%   3GPP TS 38.212 version 17.3.0 Release 17 ETSI TS 138 212 V17.3.0 (2022-09)
%   5.3.1.2 Polar encoding
A = length(a);
C = 1;
G = E;

%%% Use CA-polar
% The CRC polynomial used with CA-polar in 3GPP PUCCH channel is
% D^11 + D^10 + D^9 + D^5 + 1
crc_polynomial_pattern = [1 1 1 0 0 0 1 0 0 0 0 1];

% Determine the number of information and CRC bits.
P = length(crc_polynomial_pattern)-1;
K = ceil(A/C)+P;
E_r = floor(G/C);

% Determine the number of bits used at the input and output of the polar
% encoder kernal.
N = get_3GPP_N(K,E_r,10); % n_max = 10 is used in PUCCH channels

% Get a rate matching pattern.
[rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K,N,E_r);

% Get a sequence pattern.
Q_N = get_3GPP_sequence_pattern(N);

% Get the channel interleaving pattern
channel_interleaver_pattern = get_3GPP_channel_interleaver_pattern(E_r);

% Get an information bit pattern.
info_bit_pattern = get_3GPP_info_bit_pattern(K, Q_N, rate_matching_pattern, mode);

%%% Perform polar encoding.
%   e = CA_polar_encoder(a,crc_polynomial_pattern,info_bit_pattern,rate_matching_pattern);

% Generate the CRC bits and append them to the information bits.
G_P = get_crc_generator_matrix(A,crc_polynomial_pattern);
a_crc = mod(a*G_P,2);
b = [a, a_crc];

% Position the information and CRC bits within the input to the polar
% encoder kernal.
u = zeros(1,N);
u(info_bit_pattern) = b;

% Perform the polar encoder kernal operation.
G_N = get_G_N(N);
d = mod(u*G_N,2);

% Extract the encoded bits from the output of the polar encoder kernal.
e = d(rate_matching_pattern);

% Perform channel interleaving.
f = e(channel_interleaver_pattern);

% Testvectors
if not(exist(tc))
  mkdir(tc)
end
printf("Polar encoding: A:%d, G:%d, C:%d, P:%d, K:%d, E_r:%d, N:%d, mode:%s\n", A, G, C, P, K, E_r, N, mode);
tvwrite([tc "/" "params.txt"], [A, G, C, P, K, E_r, N]);
tvwrite([tc "/" "info_bits.txt"], a);
tvwrite([tc "/" "rate_matching_pattern.txt"], rate_matching_pattern-1);
tvwrite([tc "/" "info_bit_pattern.txt"], info_bit_pattern);
tvwrite([tc "/" "crc_polynomial_pattern.txt"], crc_polynomial_pattern);
tvwrite([tc "/" "polar_sequence.txt"], Q_N-1);
tvwrite([tc "/" "channel_interleaver_pattern.txt"], channel_interleaver_pattern-1);
tvwrite([tc "/" "info_crc_bits.txt"], b);
tvwrite([tc "/" "info_frozen_bits.txt"], u);
tvwrite([tc "/" "enc_bits.txt"], d);
tvwrite([tc "/" "rm_bits.txt"], e);
tvwrite([tc "/" "interl_bits.txt"], f);

status = 0;

