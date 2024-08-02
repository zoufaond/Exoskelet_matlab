function force = muscle_force(act, lmt, fmax, lceopt, lslack)
    % thelen 2003 passive and active length properties
    lm = lmt - lslack;
    f_gauss = 0.25;
    kpe = 5;
    epsm0 = 0.6;
    fpe = (exp(kpe*(lm / lceopt - 1)/epsm0)-1)/(exp(kpe)-1);
    flce = (exp(-(lm / lceopt - 1).^2 / f_gauss));
    force = (flce * act +  fpe) * fmax;
end