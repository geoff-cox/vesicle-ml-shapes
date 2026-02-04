function lit = sql_blob_lit(u8)
%SQL_BLOB_LIT Create SQLite BLOB literal X'...'
% u8 must be a uint8 vector.

    if isempty(u8)
        % Empty blob literal is X'' (valid)
        lit = "X''";
        return;
    end
    if ~isa(u8,'uint8')
        error('sql_blob_lit:type','Expected uint8.');
    end
    hex = upper(reshape(dec2hex(u8,2).', 1, [])); % 2 hex chars per byte
    lit = "X'" + string(hex) + "'";
end
