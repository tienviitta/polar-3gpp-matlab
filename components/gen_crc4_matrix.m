% Message
A = 12;
msg = round(rand(A, 1))
%msg = [1; 0; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1]
l_msg = length(msg); % Note! Should be divisible with step!

% The CRC polynomial
% D^4 + D^3 + 1
crc4 = [1 1 0 0 1];
l_crc = length(crc4);
poly = transpose(crc4);

% Add zeros to the end to make size step
step = 4;
poly = [poly; zeros(step-l_crc,1)];

% Add identity matrix as in
% https://en.wikipedia.org/wiki/Linear-feedback_shift_register#Matrix_forms
crc_m = zeros(step, step);
crc_m(1:end-1, 2:end) = eye(step-1);

% Add polynomial coefficients. The first one (coefficient of x^P) is not included.
% It is assumed to be always 1
crc_m(:, 1) = poly(2:end);

% Calculate matrix^step. Then the result matrix will do steps of CRC with
% one matrix-vector multiplication
crc_m = mod(crc_m^step,2);

% Compute CRC
crc = zeros(step, 1);
for i = 0:(l_msg / step) - 1
  p_msg = msg(i*step+1:(i+1)*step);
  crc = mod(crc_m * mod(crc + p_msg, 2), 2);
end
crc

% Attach CRC
r_msg = [msg; crc];
l_msg = length(r_msg);

%r_msg(7) = not(r_msg(7)); % Introduce an error

% Compute CRC
crc = zeros(step, 1);
for i = 0:(l_msg / step) - 1
  p_msg = r_msg(i*step+1:(i+1)*step);
  crc = mod(crc_m * mod(crc + p_msg, 2), 2);
end
crc

