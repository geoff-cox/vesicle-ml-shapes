function hex = hash_sha256(str)
%HASH_SHA256 SHA-256 hex digest of UTF-8 string.
    if isempty(str), str = ""; end
    bytes = uint8(unicode2native(char(string(str)), 'UTF-8'));
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(bytes);
    digest = typecast(md.digest(), 'uint8');
    hex = lower(string(reshape(dec2hex(digest)', 1, [])));
end
