function sol2 = coarsen_mesh(sol1, ratio)
% Simple resample to a coarser grid to regularize Jacobian.
    if nargin<2, ratio = 0.5; end
    x = sol1.x;  m = max(6, ceil(numel(x)*ratio));
    x2 = linspace(x(1), x(end), m);
    y2 = deval(sol1, x2);
    sol2 = bvpinit(x2, y2, sol1.parameters);
end