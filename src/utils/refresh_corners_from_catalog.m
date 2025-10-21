function [C, anyUnknown] = refresh_corners_from_catalog(C, T, MP, tol)
    % For each corner, check if a matching param exists in catalog T.
    % If yes, fill label/energy/pressure and mark solved.
    % A simple exact match on (H0_1,H0_2) is used; you can add tolerances if needed.
    anyUnknown = false;
    for i=1:4
        H1 = C.corners(i,1); H2 = C.corners(i,2);
        [solved,label,E,P] = lookup_in_catalog(T,H1,H2,MP,tol);
        C.cornerSolved(i)=solved; C.cornerLabel(i)=label;
        C.cornerEnergy(i)=E; C.cornerPressure(i)=P;
        if ~solved, anyUnknown=true; end
    end
end