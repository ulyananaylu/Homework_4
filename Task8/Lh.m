function tmp = Lh(x, t, Uiii, Uii, Ui, h)
    tmp = a(x, t) * (Uiii - 2 * Uii + Ui) / (h ^ 2) + b(x, t) * (Uiii - Ui) / (2 * h) + c(x, t) * Uii;
end
