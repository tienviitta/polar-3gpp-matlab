function polar_enc_dl(tc = 0, seed = 0)

addpath 'components'

% Seed the random number generator
rng(seed);

% Params
switch tc
  case 0
    A = 12;
    E = 48;
    EsN0 = 15;
    tc = "dl/tv0";
  case 1
    A = 65;
    E = 184;
    EsN0 = 7;
    tc = "dl/tv1";
  case 2
    % mode:repetition
    A = 134;
    E = 268;
    EsN0 = 1;
    tc = "dl/tv2";
end

% N0
N0 = 1/(10^(EsN0/10));

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

% QPSK modulation
f2 = [f,zeros(1,mod(-length(f),2))];
tx = sqrt(1/2)*(2*f2(1:2:end)-1)+1i*sqrt(1/2)*(2*f2(2:2:end)-1);

% Simulate transmission
rx = tx + sqrt(N0/2)*(randn(size(tx))+1i*randn(size(tx)));

% QPSK demodulation
f2_tilde = zeros(size(f2));
f2_tilde(1:2:end) = -4*sqrt(1/2)*real(rx)/N0;
f2_tilde(2:2:end) = -4*sqrt(1/2)*imag(rx)/N0;
f_tilde = f2_tilde(1:length(f));

% Perform polar decoding
L = 4;
min_sum = true;
a_hat = PDCCH_decoder_ref(f_tilde, A, L, min_sum, RNTI);

% Check results
errs = sum(abs(a_hat - a));
if errs
  printf("  FAIL\n");
else
  printf("  PASS\n");
end

end

