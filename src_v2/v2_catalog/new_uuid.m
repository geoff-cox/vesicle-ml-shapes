function id = new_uuid()
%NEW_UUID UUID string.
    id = string(char(java.util.UUID.randomUUID()));
end
