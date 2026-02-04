function b = bytes_per_dtype(dtype)
    switch dtype
        case "float32", b = 4;
        case "float64", b = 8;
        case "int32",   b = 4;
        case "int64",   b = 8;
        otherwise, error('bytes_per_dtype:bad','Unknown dtype.');
    end
end