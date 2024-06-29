function force = muscle_force(act, lmt, fmax, lceopt, lslack)
    lm = lmt - lslack;
    f_gauss = 0.25;
    force = (((lm / lceopt)^3) * exp(8 * lm / lceopt - 12.9) + (exp(-(lm / lceopt - 1)^2 / f_gauss)) * act) * fmax;
end