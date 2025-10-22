function [C1,C2,C3,C4] = subdivideCell(C)
    % Split into 4 children (SW, SE, NE, NW), inherit depth+1.
    am = 0.5*(C.a1 + C.a2);
    bm = 0.5*(C.b1 + C.b2);
    d1 = C.depth + 1;

    C1 = makeCell(C.a1, am, C.b1, bm, d1); % SW
    C2 = makeCell(am,  C.a2, C.b1, bm, d1); % SE
    C3 = makeCell(am,  C.a2, bm,  C.b2, d1); % NE
    C4 = makeCell(C.a1, am, bm,  C.b2, d1); % NW
end