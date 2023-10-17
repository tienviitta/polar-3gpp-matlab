function polar_enc_dl(tc = 0, seed = 0)

addpath 'components'

% Seed the random number generator
rng(seed);

% Params
switch tc
  case 0
    A = 12;
    E = 48;
    tc = "dl/tv0";
  case 1
    A = 65;
    E = 184;
    tc = "dl/tv1";
  case 2
    % mode:repetition
    A = 134;
    E = 267;
    tc = "dl/tv2";
end

% Check if testcase directory exists
if not(exist(tc))
  mkdir(tc);
end

% Generate a random block of bits
a = round(rand(1,A));

% RNTI
RNTI = randi(2, 1, 16) - 1;

% Polar Distributed-CRC-Aided encoding; Reference model
f = PDCCH_encoder_ref(a, E, RNTI, tc);

end

