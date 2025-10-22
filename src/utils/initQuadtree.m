function QT = initQuadtree(grid)
    % grid.A_lim = [amin, amax]; grid.B_lim = [bmin, bmax];
    C = makeCell(grid.A_lim(1), grid.A_lim(2), ...
                 grid.B_lim(1), grid.B_lim(2), 0);
    QT = struct('cells',[],'queue',{{C}});
end
