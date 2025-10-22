function [C1,C2,C3,C4] = subdivideCell(C)
% Split into 4 children (SW, SE, NE, NW), inherit depth+1.
    a1 = C.corners(1,1);
    a2 = C.corners(2,1);
    b1 = C.corners(2,2);
    b2 = C.corners(3,2);
    am = 0.5*(a1 + a2);
    bm = 0.5*(b1 + b2);
    d1 = C.depth + 1;

    C1 = makeCell([a1, am; b1, bm], d1);  % SW
    C2 = makeCell([am, a2; b1, bm], d1);  % SE
    C3 = makeCell([am, a2; bm, b2], d1);  % NE
    C4 = makeCell([a1, am; bm, b2], d1);  % NW
end