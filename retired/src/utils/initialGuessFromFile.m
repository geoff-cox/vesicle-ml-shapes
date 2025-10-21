function initSol = initialGuessFromFile(S, ~)
    p1 = fullfile(S.paths.mats, S.initialGuess);
    if exist(p1,'file')~=2
        p1 = fullfile('InitialShapes', S.initialGuess);
    end
    tmp = load(p1);
    initSol = tmp.Version(1).Solution;
end