function n = file_bytes(filepath)
%FILE_BYTES Return file size in bytes.
%
%   n = file_bytes(filepath)

    info = dir(filepath);
    if isempty(info) || ~isfield(info,'bytes')
        error('file_bytes:notfound','File not found: %s', filepath);
    end
    n = double(info.bytes);
end
