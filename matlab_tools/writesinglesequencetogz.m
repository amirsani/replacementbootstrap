function writesinglesequencetogz(data,file_name)

    fid = fopen(file_name,'wb');
    fwrite(fid, data, 'int');    
    fclose(fid);
    system(['gzip ' file_name]);
end
