clc
clear all

% CRC computation and attachment
% A: 12, L: 6, E: 28
% inp: [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, ]
% crc: [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1]
% ref: [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1]
%a = [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]
a = [1; 0; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1]

% The CRC polynomial used with PCCA-polar in 3GPP PUCCH channel is
% D^6 + D^5 + 1
crc6 = [1 1 0 0 0 0 1];
l_crc = length(crc6)
poly = transpose(crc6);

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
crc_m = mod(crc_m^step,2)

% Compute CRC
%crc = mod([fliplr(a), zeros(1, 32-12)] * crc_m, 2)
%crc = fliplr(mod([fliplr(a), zeros(1, 32-12)] * crc_m(:, 1:l_crc-1), 2))
%crc = mod(crc_m * [zeros(32-12, 1); a],2)
crc = mod(crc_m(1:l_crc-1, :) * [zeros(32-12, 1); a],2)

%matrix_values = uint32(zeros(step,1));
%for i = 1:step
%    row = crc_m(i,:);
%    value = 0;
%    for j=1:step
%        value = value + row(j)*2^(j-1); % Leftmost bit is interpreted as LSB
%    end
%    matrix_values(i) = value;
%end
%matrix_values

