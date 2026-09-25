function vectorEPP = enforcePositivityPatttern(vector,func)

    vectorEPP = vector;
    r=0;

    [~,signPattern] = func(zeros(size(vectorEPP)),r);
    setPos = find(signPattern == 1);
    vectorEPP(setPos) = abs(vector(setPos));
    setNeg = find(signPattern == -1);
    vectorEPP(setNeg) = -abs(vector(setNeg));

end