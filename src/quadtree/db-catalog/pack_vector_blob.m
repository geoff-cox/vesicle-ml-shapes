function [u8, dtype, n] = pack_vector_blob(vec, dtype)
%PACK_VECTOR_BLOB Pack numeric vector into uint8 bytes (little-endian).
%
% dtype: 'float32'|'float64'|'int32'|'int64'
% Returns:
%   u8   uint8 row vector
%   dtype normalized dtype string
%   n    number of elements

    arguments
        vec
        dtype (1,1) string {mustBeMember(dtype, ["float32","float64","int32","int64"])}
    end

    if isempty(vec)
        u8 = uint8([]);
        n = 0;
        return;
    end

    vec = vec(:); % column
    n = numel(vec);

    switch dtype
        case "float32"
            v = single(vec);
            u = typecast(v, 'uint8');
        case "float64"
            v = double(vec);
            u = typecast(v, 'uint8');
        case "int32"
            v = int32(vec);
            u = typecast(v, 'uint8');
        case "int64"
            v = int64(vec);
            u = typecast(v, 'uint8');
    end

    % Force little-endian representation
    if computer('endian') == 'B'
        % swap bytes within each element
        bytesPer = bytes_per_dtype(dtype);
        u = reshape(u, bytesPer, []);
        u = flipud(u);
        u = u(:);
    end

    u8 = uint8(u(:).'); % row
end