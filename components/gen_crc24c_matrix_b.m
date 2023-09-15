clc

% Random Message
%A = 1234;
%msg = round(rand(A, 1));

% Fixed message
msg = [
  1; 0; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 0; 0; 0;
  1; 0; 1; 0; 1; 1; 1; 1; 0; 0; 1; 0; 1; 1; 0; 1;
  0; 1; 1; 1; 0; 1; 1; 0; 0; 1; 1; 1; 1; 0; 1; 1;
  0; 0; 1; 0; 0; 0; 1; 1; 0; 0; 1; 1; 1; 0; 0; 1;
  1; 1; 0; 1; 1; 1; 0; 1; 0; 1; 0; 1; 1; 1; 1; 1;
  1; 0; 1; 1; 1; 1; 0; 1; 1; 1; 0; 1; 1; 0; 0; 0;
];
l_msg = length(msg);

% The CRC polynomial
crc24c = [1 1 1 0 1 0 0 0 1 0 0 0 1 1 0 1 0 1 0 0 1 1 0 1 1];
l_crc = length(crc24c);
poly = transpose(crc24c);

% Reference CRC
G_P = get_crc_generator_matrix(l_msg, crc24c);
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
crc_m = mod(crc_m ^ step,2);
%disp("enc: crc_m:")
%disp(crc_m)

% Compute CRC in 'step' size blocks
n_pad = ceil(l_msg / step) * step - l_msg;
printf("enc: pad: %d\n", n_pad);
z_msg = [zeros(n_pad, 1); msg];
l_msg = length(z_msg);
crc_s = zeros(step, 1);
for i = 0:(l_msg / step) - 1
  p_msg = z_msg(i*step+1:(i+1)*step);
  printf("enc: p_msg(%d): ", i);
  printf("%d", fliplr(transpose(p_msg)));
  printf("\n");
  crc_p = xor(crc_s, p_msg);
  printf("enc: crc_p(%d): ", i);
  printf("%d", fliplr(transpose(crc_p)));
  printf("\n");
  %crc_s = mod(crc_m * crc_p, 2); % Note! Zero rows?!
  for j = 1:step
    crc_v = and(transpose(crc_m(j,:)), crc_p);
    crc_s(j) = mod(sum(crc_v), 2);
%    printf("enc: crc_m(%d,:): ", j);
%    printf("%d", fliplr(crc_m(j,:)));
%    printf("\n");
%    printf("enc: crc_v(%d,:): ", j);
%    printf("%d", flipud(crc_v));
%    printf("\n");
  endfor
  printf("enc: crc_s(%d): ", i);
  printf("%d", fliplr(transpose(crc_s)));
  printf("\n");
endfor
printf("enc: crc_s: ");
printf("%d", fliplr(transpose(crc_s(1:l_crc-1))));
printf("\n");

% Attach CRC
r_msg = [msg; crc_s(1:l_crc-1)];
l_msg = length(r_msg);

%i_err = randi(A);
%r_msg(i_err) = not(r_msg(i_err)); % Note! Introduce a bit error to the msg!

% Compute CRC in 'step' size blocks
n_pad = ceil(l_msg / step) * step - l_msg
printf("dec: pad: %d\n", n_pad);
r_msg = [zeros(n_pad, 1); r_msg];
l_msg = length(r_msg);
crc_s = zeros(step, 1);
for i = 0:(l_msg / step) - 1
  p_msg = r_msg(i*step+1:(i+1)*step);
  printf("dec: p_msg(%d): ", i);
  printf("%d", fliplr(transpose(p_msg)));
  printf("\n");
  crc_p = xor(crc_s, p_msg);
  printf("dec: crc_p(%d): ", i);
  printf("%d", fliplr(transpose(crc_p)));
  printf("\n");
  %crc_s = mod(crc_m * crc_p, 2); % Note! Zero rows?!
  for j = 1:step
%    crc_s(j) = mod(sum(and(transpose(crc_m(j,:)), crc_p)), 2);
    crc_v = and(transpose(crc_m(j,:)), crc_p);
    crc_s(j) = mod(sum(crc_v), 2);
  endfor
  printf("dec: crc_s(%d): ", i);
  printf("%d", fliplr(transpose(crc_s)));
  printf("\n");
endfor
printf("dec: crc_s: ");
printf("%d", fliplr(transpose(crc_s(1:l_crc-1))));
printf("\n");
if sum(crc_s(1:l_crc-1)) == 0
  printf("dec: PASS\n");
else
  printf("dec: FAIL\n");
endif

% LUT for DSP (bit-by-bit) PUCCH encoder implementation
matrix_values = uint32(zeros(step,1));
for i = 1:step
    row = crc_m(i,:);
    value = 0;
    for j=1:step
        value = value + row(j)*2^(j-1); % Note! Leftmost bit is interpreted as LSB!
    end
    matrix_values(i) = value;
end
disp("lut:")
disp(dec2hex(matrix_values))

