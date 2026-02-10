function s = iso_utc(dt)
%ISO_UTC Convert datetime to ISO-8601 Z string.
    if isempty(dt)
        s = "";
        return;
    end
    if isempty(dt.TimeZone)
        dt.TimeZone = 'UTC';
    else
        dt = datetime(dt, 'TimeZone','UTC');
    end
    s = string(datestr(dt, 'yyyy-mm-ddTHH:MM:SS')) + "Z";
end
