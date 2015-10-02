function [x] = read_sequences(file_name, n, S)

x = zeros(S,n);

% Read all the sequences
system(['cp ' file_name ' temp.dat.gz']);
system('gunzip temp.dat.gz');
fid = fopen('temp.dat');

for t=1:n
    x(:,t) = fread(fid, S, 'int');
end

system('rm temp.dat');
fclose(fid);

end
