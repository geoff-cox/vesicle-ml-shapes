function tf = is_hash_in_catalog(T, hash)
    if isempty(T); tf = false; return; end
    tf = any(T.hash == string(hash));
end