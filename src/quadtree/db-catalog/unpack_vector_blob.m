function vec = unpack_vector_blob(u8, dtype, n)
%UNPACK_VECTOR_BLOB Decode uint8 bytes (little-endian) into numeric vector.

    arguments
        u8 (1,:) uint8
        dtype (1,1) string {mustBeMember(dtype, ["float32","float64","int32","int64"])}
        n (1,1) double {mustBeInteger, mustBeNonnegative}
    end

    if n == 0
        vec = zeros(0,1);
        return;
    end

    bytesPer = bytes_per_dtype(dtype);
    if numel(u8) ~= n * bytesPer
        error('unpack_vector_blob:size','Byte length mismatch.');
    end

    u = u8(:);

    % Convert from little-endian into native if needed
    if computer('endian') == 'B'
        u = reshape(u, bytesPer, []);
        u = flipud(u);
        u = u(:);
    end

    switch dtype
        case "float32"
            vec = typecast(u, 'single');
        case "float64"
            vec = typecast(u, 'double');
        case "int32"
            vec = typecast(u, 'int32');
        case "int64"
            vec = typecast(u, 'int64');
    end

    vec = vec(:);
end