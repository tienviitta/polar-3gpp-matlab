clear all
close all
clc

addpath 'components'

% Read message
msg = dlmread("dl/tv0/info_bits.txt");

l_msg = length(msg);

% The CRC polynomial used in 3GPP PBCH and PDCCH channel is
% D^24 + D^23 + D^21 + D^20 + D^17 + D^15 + D^13 + D^12 + D^8 + D^4 + D^2 + D + 1
crc_polynomial_pattern = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
l_crc = length(crc_polynomial_pattern);
poly = transpose(crc_polynomial_pattern);

% Reference CRC
G_P = get_crc_generator_matrix(l_msg, crc_polynomial_pattern);
crc_ref = mod(transpose(msg) * G_P, 2);
printf("ref: crc_s: ");
printf("%d", transpose(fliplr(crc_ref)));
printf("\n");

% Add zeros to the end to make size step
step = 32;
poly = [poly; zeros(step-l_crc,1)];

% Add identity matrix as in
% https://en.wikipedia.org/wiki/Linear-feedback_shift_register#Matrix_forms
crc_m = zeros(step, step);
crc_m(1:end-1, 2:end) = eye(step-1);

% Add polynomial coefficients. The first one (coefficient of x^P) is not included.
% It is assumed to be always 1
crc_m(1:step-1,1) = poly(2:end);

% Calculate matrix^step. Then the result matrix will do steps of CRC with
% one matrix-vector multiplication
crc_b = mod(crc_m ^ step,2);
step_rem = rem(l_msg, step);
crc_r = mod(crc_m ^ step_rem,2);
%crc_r(:,end-(32-step_rem-1):end) = 0;

% Compute CRC in 'step' size blocks
crc_s = zeros(step, 1);
for i = 0:floor(l_msg / step) - 1
  p_msg = msg(i*step+1:(i+1)*step);
  printf("enc: p_msg(%d): ", i);
  printf("%d", fliplr(transpose(p_msg)));
  printf("\n");
  crc_p = xor(crc_s, p_msg);
  printf("enc: crc_p(%d): ", i);
  printf("%d", fliplr(transpose(crc_p)));
  printf("\n");
  for j = 1:step
    crc_v = and(transpose(crc_b(j,:)), crc_p);
    crc_s(j) = mod(sum(crc_v), 2);
  endfor
  printf("enc: crc_s(%d): ", i);
  printf("%d", fliplr(transpose(crc_s)));
  printf("\n");
endfor

% Compute remainder
p_msg = [msg(floor(l_msg/step)*step+1:floor(l_msg/step)*step+rem(l_msg, step)); zeros(32-rem(l_msg,32), 1)];
printf("enc: p_msg(%d): ", floor(l_msg/step));
printf("%d", fliplr(transpose(p_msg)));
printf("\n");
crc_p = xor(crc_s, p_msg);
printf("enc: crc_p(%d): ", floor(l_msg/step));
printf("%d", fliplr(transpose(crc_p)));
printf("\n");
for j = 1:step
  crc_v = and(transpose(crc_r(j,:)), crc_p);
  crc_s(j) = mod(sum(crc_v), 2);
endfor
printf("enc: crc_s(%d): ", floor(l_msg/step));
printf("%d", fliplr(transpose(crc_s)));
printf("\n");
printf("enc: crc_s: ");
printf("%d", fliplr(transpose(crc_s(1:l_crc-1))));
printf("\n");

if all(crc_ref == transpose(crc_s(1:l_crc-1)))
  printf("enc: PASS\n");
else
  printf("enc: FAIL\n");
endif

%% Attach CRC
%r_msg = [msg; crc_s(1:l_crc-1)];
%l_msg = length(r_msg);
%
%%i_err = randi(A);
%%r_msg(i_err) = not(r_msg(i_err)); % Note! Introduce a bit error to the msg!
%
%% Compute CRC in 'step' size blocks
%n_pad = ceil(l_msg / step) * step - l_msg
%printf("dec: pad: %d\n", n_pad);
%r_msg = [zeros(n_pad, 1); r_msg];
%l_msg = length(r_msg);
%crc_s = zeros(step, 1);
%for i = 0:(l_msg / step) - 1
%  p_msg = r_msg(i*step+1:(i+1)*step);
%  printf("dec: p_msg(%d): ", i);
%  printf("%d", fliplr(transpose(p_msg)));
%  printf("\n");
%  crc_p = xor(crc_s, p_msg);
%  printf("dec: crc_p(%d): ", i);
%  printf("%d", fliplr(transpose(crc_p)));
%  printf("\n");
%  %crc_s = mod(crc_b * crc_p, 2); % Note! Zero rows?!
%  for j = 1:step
%%    crc_s(j) = mod(sum(and(transpose(crc_b(j,:)), crc_p)), 2);
%    crc_v = and(transpose(crc_b(j,:)), crc_p);
%    crc_s(j) = mod(sum(crc_v), 2);
%  endfor
%  printf("dec: crc_s(%d): ", i);
%  printf("%d", fliplr(transpose(crc_s)));
%  printf("\n");
%endfor
%printf("dec: crc_s: ");
%printf("%d", fliplr(transpose(crc_s(1:l_crc-1))));
%printf("\n");

%if sum(crc_s(1:l_crc-1)) == 0
%  printf("dec: PASS\n");
%else
%  printf("dec: FAIL\n");
%endif

% LUT for DSP (bit-by-bit) PUCCH encoder implementation
matrix_values = uint32(zeros(step,1));
for i = 1:step
    row = crc_b(i,:);
    value = 0;
    for j=1:step
        value = value + row(j)*2^(j-1); % Note! Leftmost bit is interpreted as LSB!
    end
    matrix_values(i) = value;
end
disp("lut:")
disp(dec2hex(matrix_values))

mask = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0];
matrix_values = uint32(zeros(step,1));
for i = 1:step
    row = crc_r(i,:);
    %row = and(crc_r(i,:), mask);
    value = 0;
    for j=1:step
        value = value + row(j)*2^(j-1); % Note! Leftmost bit is interpreted as LSB!
    end
    matrix_values(i) = value;
end
disp("lut:")
disp(dec2hex(matrix_values))

