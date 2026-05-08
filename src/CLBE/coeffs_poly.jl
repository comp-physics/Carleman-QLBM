function coeffs_all(C, f)
        dum = sympy.Poly(C, f)
        dum_coeff = dum.coeffs()
	vars = [prod(dum.gens.^I) for I in Iterators.product((0:d for d in dum.degree.(dum.gens))...)]
        vars = reshape(vars,(length(vars), 1))
	return vars
end

function coeffs(C, f)
        dum = sympy.Poly(C, f)
	 ms = dum.coeffs()
	 vs = [prod(dum.gens.^i) for i in dum.monoms()]
	return vs, ms
end
