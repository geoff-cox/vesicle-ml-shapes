function hex = file_sha256(filepath, chunkBytes)
%FILE_SHA256 Compute SHA-256 hex digest of a file.
%
%   hex = file_sha256(filepath)
%   hex = file_sha256(filepath, chunkBytes)
%
% Uses streaming read so it works on large HDF5 files.

    if nargin < 2 || isempty(chunkBytes)
        chunkBytes = 16 * 1024 * 1024; % 16 MB
    end
    if chunkBytes < 1024
        error('file_sha256:chunk','chunkBytes too small.');
    end

    fid = fopen(filepath, 'rb');
    if fid < 0
        error('file_sha256:open','Could not open file: %s', filepath);
    end

    cleaner = onCleanup(@() fclose(fid));

    md = java.security.MessageDigest.getInstance('SHA-256');

    while true
        data = fread(fid, chunkBytes, '*uint8');
        if isempty(data)
            break;
        end
        md.update(data);
    end

    digest = typecast(md.digest(), 'uint8');
    hex = lower(string(reshape(dec2hex(digest)', 1, [])));
end
