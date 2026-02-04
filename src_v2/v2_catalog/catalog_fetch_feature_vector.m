function vec = catalog_fetch_feature_vector(db, solution_id, pack_id, vector_name)
%CATALOG_FETCH_FEATURE_VECTOR Fetch and decode a vector from feature_vectors.

    sql = "SELECT dtype,length,blob FROM feature_vectors WHERE " + ...
          "solution_id=" + sql_lit(solution_id) + " AND pack_id=" + sql_lit(pack_id) + ...
          " AND vector_name=" + sql_lit(vector_name) + " LIMIT 1;";

    rows = fetch(db, sql);
    if isempty(rows)
        vec = [];
        return;
    end

    dtype = string(rows{1,1});
    n = double(rows{1,2});
    blob = rows{1,3};

    if isa(blob,'uint8')
        u8 = blob(:).';
    elseif ischar(blob) || isstring(blob)
        % If it comes back as hex text for some reason, try to parse:
        u8 = hex_to_uint8(string(blob));
    else
        error('catalog_fetch_feature_vector:blobtype','Unexpected blob type: %s', class(blob));
    end

    vec = unpack_vector_blob(uint8(u8), dtype, n);
end

function u8 = hex_to_uint8(hex)
    hex = regexprep(char(hex), '\s+', '');
    if startsWith(hex,"X'") || startsWith(hex,"x'")
        hex = hex(3:end-1);
    end
    if mod(numel(hex),2)~=0
        error('hex_to_uint8:bad','Hex string length must be even.');
    end
    u8 = uint8(sscanf(hex,'%2x').');
end
