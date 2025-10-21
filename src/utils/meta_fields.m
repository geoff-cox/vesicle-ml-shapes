function [label, E, P] = meta_fields(m)
    label = "";
    E = NaN; P = NaN;
    if ~isstruct(m), return; end
    if isfield(m,'label'); label = string(m.label); end
    if isfield(m,'E');     E     = double(m.E);    end
    if isfield(m,'P');     P     = double(m.P);    end
    % Fallbacks (older names)
    if isnan(E) && isfield(m,'energy'),   E = double(m.energy);   end
    if isnan(P) && isfield(m,'pressure'), P = double(m.pressure); end
    if label=="" && isfield(m,'shapeLabel'), label = string(m.shapeLabel); end
end