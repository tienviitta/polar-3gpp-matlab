function polar_enc_vec(tc = 0, seed = 0)

%clc
%clear all
%close all

addpath 'components'

%%% Params
switch tc
  case 0
    % mode:shortening
    A = 123;
    E = 184;
    tc = "tv0";
  case 1
    % mode:puncturing
    A = 65;
    E = 184;
    tc = "tv1";
  case 2
    % mode:repetition
    A = 134;
    E = 267;
    tc = "tv2";
end
%seed = 0;

% Polar CRC-aided encoding; Reference model
polar_enc_ref(A, E, tc, seed)

end

