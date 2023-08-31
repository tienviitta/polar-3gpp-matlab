% Message
A = 3*17;
msg = round(rand(A, 1));
l_msg = length(msg);

% The CRC polynomial used with CA-polar in 3GPP PUCCH channel is
% D^11 + D^10 + D^9 + D^5 + 1
crc11 = [1 1 1 0 0 0 1 0 0 0 0 1];
l_crc = length(crc11);
poly = transpose(crc11);

% Reference CRC
G_P = get_crc_generator_matrix(A, crc11);
crc_ref = mod(transpose(msg) * G_P, 2);
disp("ref: crc_s:")
disp(crc_ref)

% Add zeros to the end to make size step
step = 16;
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
disp("enc: crc_m:")
disp(crc_m)

% Compute CRC in 'step' size blocks
n_pad = ceil(l_msg / step) * step - l_msg
z_msg = [zeros(n_pad, 1); msg];
l_msg = length(z_msg);
crc_s = zeros(step, 1);
for i = 0:(l_msg / step) - 1
  p_msg = z_msg(i*step+1:(i+1)*step);
  crc_s = mod(crc_m * mod(crc_s + p_msg, 2), 2); % Note! Zero rows?!
end
disp("enc: crc_s:")
disp(transpose(crc_s))

% Attach CRC
r_msg = [msg; crc_s(1:l_crc-1)];
l_msg = length(r_msg);

%i_err = randi(A);
%r_msg(i_err) = not(r_msg(i_err)); % Note! Introduce a bit error to the msg!

% Compute CRC in step size blocks
n_pad = ceil(l_msg / step) * step - l_msg
r_msg = [zeros(n_pad, 1); r_msg];
l_msg = length(r_msg);
crc_s = zeros(step, 1);
for i = 0:(l_msg / step) - 1
  p_msg = r_msg(i*step+1:(i+1)*step);
  crc_s = mod(crc_m * mod(crc_s + p_msg, 2), 2); % Note! Zero rows?!
end
disp("dec: crc_s:")
disp(transpose(crc_s))

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

