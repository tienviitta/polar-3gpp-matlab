crc24c = zeros(25,1);
onesList = [24 23 21 20 17 15 13 12 8 4 2 1 0];
for i = onesList
    crc24c(i+1) = 1;
end


pol = flipud(crc24c(1:25));


sz = 32;
pol = [pol; zeros(sz-25,1)];    %add zeros to the end to make size 32

crc_m = zeros(sz,sz);

% add identity matrix as in https://en.wikipedia.org/wiki/Linear-feedback_shift_register#Matrix_forms
crc_m(1:end-1, 2:end) = eye(sz-1);

% add polynomial coefficients. The first one (coefficient of x^24) is not included.
% it is assumed to be always 1
crc_m(1:sz-1,1) = pol(2:end);

% Calculate matrix^32. Then the result matrix will do 32 steps of CRC with
% one matrix-vector multiplication
crc_m = mod(crc_m^sz,2)


matrix_values = uint32(zeros(32,1));
for i = 1:32
    row = crc_m(i,:);
    value = 0;
    for j=1:32
        value = value + row(j)*2^(j-1); %Leftmost bit is interpreted as LSB
    end
    matrix_values(i) = value;
end
dec2hex(matrix_values)

