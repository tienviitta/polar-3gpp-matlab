function a_hat = PDCCH_decoder_ref(f_tilde, A, L, min_sum, RNTI, tc)
% PDCCH_DECODER Polar decoder for the Physical Downlink Control Channel (PDCCH) of 3GPP New
% Radio, as defined in Section 7.3 of TS38.212. Implements the zero-
% padding to increase the length of short payloads to 12 bits of Section 7.3.1,
% the Cyclic Redudancy Check (CRC) attachment of Section 7.3.2, the channel
% coding of Section 7.3.3 and the rate matching of Section 7.3.4. Note that
% this code does not implement the DCI bit sequence generation of Section
% 7.3.1.
%   a_hat = PDCCH_DECODER(f_tilde, A, L, min_sum) decodes the encoded LLR sequence
%   f_tilde, in order to obtain the recovered information bit sequence
%   a_hat.
%
%   f_tilde should be a real row vector comprising E number of Logarithmic
%   Likelihood Ratios (LLRS), each having a value obtained as LLR =
%   ln(P(bit=0)/P(bit=1)). The first LLR corresponds to f_0 from Section
%   7.3.4 of TS38.212, while the last LLR corresponds to
%   f_E-1.
%
%   A should be an integer scalar. It specifies the number of bits in the
%   information bit sequence, where A should be in the range 1 to 140.
%
%   L should be a scalar integer. It specifies the list size to use during
%   Successive Cancellation List (SCL) decoding.
%
%   min_sum shoular be a scalar logical. If it is true, then the SCL
%   decoding process will be completed using the min-sum approximation.
%   Otherwise, the log-sum-product will be used. The log-sum-product gives
%   better error correction capability than the min-sum, but it has higher
%   complexity.
%
%   RNTI should be a binary row vector comprising 16 bits, each having the
%   value 0 or 1. If this parameter is omitted, then ones(1,16) will be
%   used for the RNTI. The first bit corresponds to x_rnti,0 from Section
%   7.3.2 of TS38.212, while the last bit corresponds to x_rnti,15.
%
%   a_hat will be a binary row vector comprising A number of bits, each
%   having the value 0 or 1. The first output bit corresponds to a_0 from
%   Section 7.3.1 of TS38.212, while the last output bit corresponds
%   to a_A-1.
%
%   See also PDCCH_ENCODER
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

if A == 0
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no less than 1.');
end
if A > 140
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no greater than 140.');
end

if nargin == 4
    RNTI = ones(1,16);
end
if length(RNTI) ~= 16
    error('RNTI length should be 16');
end

E = length(f_tilde);
if E > 8192
    error('polar_3gpp_matlab:UnsupportedBlockLength','E should be no greater than 8192.');
end

% The CRC polynomial used in 3GPP PBCH and PDCCH channel is
% D^24 + D^23 + D^21 + D^20 + D^17 + D^15 + D^13 + D^12 + D^8 + D^4 + D^2 + D + 1
crc_polynomial_pattern = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];

% The CRC has P bits. P-min(P2,log2(L)) of these are used for error
% detection, where L is the list size. Meanwhile, min(P2,log2(L)) of
% them are used to improve error correction. So the CRC needs to be
% min(P2,log2(L)) number of bits longer than CRCs used in other codes,
% in order to achieve the same error detection capability.
P = length(crc_polynomial_pattern)-1;
P2 = 3;

% Determine the number of information and CRC bits.
if A < 12
    % a has been padded with zeros to increase its length to 12
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

%if A < 12
%    % We know that a has been padded with zeros
%    a_tilde = [NaN(1,A),zeros(1,12-A)];
%
%    % Perform Distributed-CRC-Aided polar decoding.
%    a_hat = DS1CKA_polar_decoder(f_tilde,crc_polynomial_pattern,RNTI,crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern,mode,L,min_sum,P2,a_tilde);
%
%    if ~isempty(a_hat)
%        % Remove the padding
%        a_hat = a_hat(1:A);
%    end
%else
%    % Perform Distributed-CRC-Aided polar decoding.
%    a_hat = DS1CA_polar_decoder(e_tilde, crc_polynomial_pattern, crc_scrambling_pattern, crc_interleaver_pattern, info_bit_pattern, rate_matching_pattern, mode, L, min_sum, P2)
%    a_hat = DS1CA_polar_decoder(f_tilde, crc_polynomial_pattern, RNTI,                   crc_interleaver_pattern, info_bit_pattern, rate_matching_pattern, mode, L, min_sum, P2);
%end

% This global variable is used by the minstar and phi functions.
global approx_minstar
approx_minstar=min_sum;

%% Characterise the interleaved CRC

% Get the CRC generator matrix, which has dimensions A+P by P.
G_P = get_crc_generator_matrix(A+P,crc_polynomial_pattern);

% Extend the CRC generator matrix by append an identity matrix to
% represent the CRC bits, giving dimenstions K by P.
G_P2 = [G_P(P+1:end,:);eye(P)];

% Interleave the rows of the extended CRC generator matrix, according to
% the CRC interleaver.
G_P3 = G_P2(crc_interleaver_pattern,:);

% Determine where the last 1-valued bit appears in each column. When the
% SCL decoding process reaches the corresponding interleaved CRC bit, we
% will terminate the decoding process if all list entries have failed the
% checks assocated with at least one of this or the preceeding CRC bits.
last_one_index = zeros(1,P);
for p = 1:P
    last_one_index(p) = find(G_P3(:,p) == 1, 1, 'last');
end

% Extend the scrambling pattern to match the length of the CRC
extended_crc_scrambling_pattern = [zeros(1,P-length(RNTI)), RNTI];

%% Rate matching
if strcmp(mode,'repetition')
    % LLRs for repeated bits are added together.
    d_tilde = zeros(1,N);
    for i=1:E
        d_tilde(rate_matching_pattern(i)) = d_tilde(rate_matching_pattern(i)) + f_tilde(i);
    end
else
    if strcmp(mode,'puncturing')
        % Zero valued LLRs are used for punctured bits, because the decoder
        % doesn't know if they have values of 0 or 1.
        d_tilde = zeros(1,N);
    elseif strcmp(mode,'shortening')
        % Infinite valued LLRs are used for shortened bits, because the
        % decoder knows that they have values of 0.
        d_tilde = inf(1,N);
    else
        error('Unknown rate matching mode');
    end

    d_tilde(rate_matching_pattern) = f_tilde;
end

%% Perform the SCL polar decoder kernal operation.
% This is achieved according to the algorithm described in
% http://ieeexplore.ieee.org/abstract/document/7114328/

global bits % Matrix to store bit throughout the polar code graph
global bits_updated % Matrix to identify which bits in the polar code graph have been calculated so far
global llrs % Matrix to store LLRs throughout the polar code graph
global llrs_updated % Matrix to identify which LLRs in the polar code graph have been calculated so far

bits = zeros(N, log2(N)+1); % Initialse all bits to zero. The left-most column corresponds to the decoded information, CRC and frozen bits
bits_updated = [~info_bit_pattern',false(N,log2(N))]; % The zero values that have initialised the frozen bits are known to be correct
llrs = [zeros(N,log2(N)),d_tilde']; % Initialse the LLRs. The right-most column corresponds to the received LLRS
llrs_updated = [false(N,log2(N)),true(N,1)]; % The received LLRs have been updated.

PM = zeros(1,1,1); % Initialise the path metrics
L_prime = 1; % Initialise the list size to 1. This will grow as the decoding proceeds

% We will calculate the CRC checksums alongside the SCL decoding process.
% Initialise the checksums with P number of 1s.
crc_checksums = mod(sum(G_P(1:P,:)),2)';
% We will keep track of whether any of the checks associated with the CRC
% bits have failed.
crc_okay = true;
% We need a counter to keep track of which information or CRC bit we are
% working on.
i2=1;
% We will return an empty vector if all list entries have failed the checks
% assocated with at least one of the CRC bits.
a_hat = [];

% Consider each bit in turn
for i = 1:N
    % Make recursive function calls to perform the XOR, g and f functions
    % necessary to obtain the corresponding LLR
    update_llr(i,1);

    if info_bit_pattern(i) == 0 % Frozen bit
        PM = phi(PM, llrs(i,1,:), 0);
    else % Information or CRC bit
        % Double the list size, using 0-valued bits for the first half and 1-valued bits for the other half
        PM = cat(3,phi(PM, llrs(i,1,:), 0), phi(PM, llrs(i,1,:), 1));
        llrs = cat(3,llrs,llrs);
        bits = cat(3,bits,bits);
        bits(i,1,1:L_prime) = 0;
        bits(i,1,L_prime+1:2*L_prime) = 1;
        bits_updated(i,1) = true;

        % We use the interleaved CRC generator matrix to update the CRC
        % check sums whenever an information or CRC bit adopts a value of
        % 1.
        crc_checksums = cat(3,crc_checksums,mod(crc_checksums+repmat(G_P3(i2,:)',[1 1 L_prime]),2));
        % We need to keep track of whether any of the checks associated
        % with the previous CRC bits have failed.
        crc_okay = cat(3,crc_okay,crc_okay);

        % If the list size has grown above L, then we need to find and keep only the best L entries in the list
        L_prime = size(bits,3);
        if L_prime > L
            [~,max_indices] = sort(PM,3);
            PM = PM(:,:,max_indices(1:L));
            bits = bits(:,:,max_indices(1:L));
            llrs = llrs(:,:,max_indices(1:L));
            crc_checksums = crc_checksums(:,:,max_indices(1:L));
            crc_okay = crc_okay(:,:,max_indices(1:L));
            L_prime = L;
        end

        % We check the corresponding CRC checksums whenever we reach the
        % last 1-valued bit in a column of the interleaved CRC generator
        % matrix.
        check_crc_bits = find(last_one_index == i2);
        for crc_bit_index = 1:length(check_crc_bits)
            for list_index = 1:L_prime
                % The checksum should equal the value of the corresponding
                % CRC scrambling pattern bit. If not, then the CRC
                % check has failed. Note that we should not prune these
                % entries from the list, even though we know that they will
                % fail the CRC. We should continue the decoding of this
                % list entries, otherwise we will damage the error
                % detection capability of the CRC.
                if crc_checksums(check_crc_bits(crc_bit_index),1,list_index) ~= extended_crc_scrambling_pattern(check_crc_bits(crc_bit_index))
                    % We keep track of the failing check.
                    crc_okay(1,1,list_index) = false;
                end
            end
        end

        % We terminate the decoding process if all list entries have failed
        % the checks assocated with at least one of this or the preceeding
        % CRC bits.
        if ~any(crc_okay)
            return;
        end

        % Increment the counter of information and CRC bits
        i2 = i2+1;
    end
end

%% Information bit extraction
tvwrite([tc "/" "dec_params.txt"], [A, P, K, E, N]);
tvwrite([tc "/" "dec_info_bit_pattern.txt"], info_bit_pattern);
d_tilde(d_tilde == Inf) = 255.0;
tvwrite([tc "/" "dec_llrs.txt"], d_tilde);

% We use the list entry with a passing CRC that has the best metric. But we
% only consider the best min(L,2^P2) entries in the list, to avoid
% degrading the CRC's error detection capability.
[~,max_indices] = sort(PM,3);
b_hat = zeros(1,K);
for list_index = 1:min(L,2^P2)
    % Consider the next best list entry.

    % We already checked the CRC during the SCL decoding process.
    if crc_okay(max_indices(list_index))
        u_hat = bits(:,1,max_indices(list_index))';

        % Extract the information bits from the output of the polar decoder
        % kernal.
        c_hat = u_hat(info_bit_pattern);

        % Deinterleave the information and CRC bits.
        b_hat(crc_interleaver_pattern) = c_hat;

        % Remove the CRC and output the information bits.
        a_hat = b_hat(1:end-P);
        return;
    end
end

