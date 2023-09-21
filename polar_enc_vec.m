clc
clear all
close all

%%% Params
A = 85;
E = 123;
seed = 0;

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
b = [a, mod(a*G_P,2)];

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

