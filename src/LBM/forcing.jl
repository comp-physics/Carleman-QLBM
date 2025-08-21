using DelimitedFiles
using LinearAlgebra
import Random
#
include("lbm_cons.jl")

Random.seed!(1234)

cs2 = lbm_param()
a, b, c, d, D, Q, w, e = lbm_cons()
#
QCFD_HOME = ENV["QCFD_HOME"] * "src/"
#---read the prepared random wave vector
function read_k()
#    k_dat = readdlm("../f90/k.dat")
    k_dat = readdlm(QCFD_HOME * "LBM/f90/k.dat")
    nk = k_dat[1,1] # number of "kkx" ("kky" or "kkz")
    kav = k_dat[1,2] # forcing factor
    #
    k_dat = replace(k_dat, ""=>NaN)
    k_dat = Vector{Float64}(vec(k_dat[2:end,:]')) # Note that Julia array indexing is column based.
    kkx = k_dat[1:nk]
    kky = k_dat[nk+2:2*nk+1]
    kkz = k_dat[2*nk+3:end-1]
    #
    return nk, kav, kkx, kky, kkz
end
#---.
nk, kav, kkx, kky, kkz =  read_k()
#= The external forcing is implemented following Eq.28 of "https://doi.org/10.1023/A:1010414013942".
We assume "C=0" in Eq.28. "B=forcing".
=#

#
function forcing_lbm(fx_random, fy_random, e, w)
    fm_forcing = Vector{Float64}(undef, Q)
    forcing_random = [fx_random, fy_random] 
    if size(e, 2) < 2
        e_temp = e
        e = zeros(Q, 2)
        e[:,1] = e_temp
    end
    for i=1:Q
	fm_forcing[i] = w[i]/cs2* dot(forcing_random, e[i, :])*dt 
    end
    #
    return fm_forcing
end
#
function forcing_random(nx, ny, nz, force_factor)
      coef1, coef2, kk, phase, fact = fconst_coefs_hel(force_factor, kkx, kky, kkz, nk, kav) 
#      println("kk", kk)
      fx_domain = zeros(nx, ny) 
      fx_domain .= NaN 
      fy_domain = zeros(nx, ny) 
      fy_domain .= NaN 
      for i = 1:nx
	for j = 1:ny
	  forcing_random_dum = forcing_hel(2*pi*i/nx, 2*pi*j/ny, 0, coef1, coef2, kk, phase, fact)[1:2] 
	  fx_domain[i,j] = forcing_random_dum[1]
	  fy_domain[i,j] = forcing_random_dum[2]
        end
      end
      #
      return fx_domain, fy_domain
end

#function random_forcing()
function forcing_hel(x, y, z, coef1, coef2, kk, phase, fact)
    # adapted from "forcing_hel", l.1854 of https://github.com/pencil-code/pencil-code/blob/master/src/forcing.f90 
    # See Eq.226 of http://pencil-code.nordita.org/doc/manual.pdf. The dimension of f(x,t), m/s^2 is given by the normalization factor N=f_/0 c_s (k c_s /delta t)^{1/2} 
    #  Add helical forcing function, using a set of precomputed wavevectors.
    #  The relative helicity of the forcing function is determined by the factor
    #  sigma, called here also relhel. If it is +1 or -1, the forcing is a fully
    #  helical Beltrami wave of positive or negative helicity. For |relhel| < 1
    #  the helicity less than maximum. For relhel=0 the forcing is nonhelical.
    #  The forcing function is now normalized to unity (also for |relhel| < 1).
    #
    #---calculate fx, fy, fz
#    fx, fy, fz, coef1, coef2 = forcing_coefs_hel(x, y, z, force_factor)
 
    fx, fy, fz = fxyz_coefs_hel(x, y, z, kk, phase, fact)
#    fx, fy, fz = forcing_coefs_hel(x, y, z, force_factor)
    #
    #
   # profyz=profy_ampl(m)*profz_ampl(n)
   # fxyz=fx(l1:l2)*fy(m)*fz(n)
   fxyz = fx .* fy .* fz # N*exp{} part of Eq.226: exp(a+b) = exp(a) * exp(b)  
    forcing_rhs = Vector{Any}(undef, 3)
    for j = 1:3
#	forcing_rhs[:,j] = force_ampl*fda[j]*cos(omega_ff*t) &
#	*real(cmplx(coef1[j],profx_hel*profyz_hel_coef2[j])*fxyz) # l.2009 of "forcing.f90"
#	forcing_rhs[:,j] = fda[j]*real(cmplx(coef1[j], coef2[j])*fxyz) # l.2009 of "forcing.f90"
forcing_rhs[j] = real((coef1[j] .+ coef2[j]*im) .* fxyz) # l.2006 of "forcing.f90" at version fa318f9e4a5283e7ed2ad6ffd6db2a7015e2220d. In the FORTRAN code, " force_ampl*fda(j)*cos(omega_ff*t)=1" and "profx_hel*profyz_hel_coef2(j) = 0" (coef2[j] = 0). 
    end				  
    #---Check NaN
    #
 forcing_rhs = collect(Iterators.flatten(forcing_rhs))
    #
    return forcing_rhs
end

function  fconst_coefs_hel(force_factor,kkx,kky,kkz,nk,kav) 
    #  This routine can be called with any values of kkx,kky,kkz (nk, kav, kkx, kky, kkz are from "read_k")
#  to produce coef1,coef2,coef3,kk,phase,fact and fda.
# "[1]" must be included after "kkx[ik]", "kky[ik]", "kkz[ik]"   
phase = pi*(2. * rand(1) .- 1.) # rand(1) generate a vector, therefore using ".-" instead of "-"
    ik= Int.(floor.(nk * (.9999 * rand(1)) .+ 1.))
    if lscale_kvector_tobox
	kx0=kkx[ik][1] .* (2. .* pi ./ Lxyz[1])
	kx = kx0
	ky=kky[ik][1] .* (2. .* pi ./ Lxyz[2])
	if Lxyz[3] != 0.
	  kz=kkz[ik][1] .* (2. .* pi ./ Lxyz[3])
	else
	    kz=kkz[ik][1]
	end
	pi_over_Lx = pi / Lxyz[1]
    else
	kx0 = kkx[ik][1]
	kx = kx0
	ky = kky[ik][1]
	kz = kkz[ik][1]
	pi_over_Lx = 0.5
    end

    #  compute k^2 and output wavenumbers
      k2 = kx .^ 2 .+ ky .^ 2 .+ kz .^ 2
      k = sqrt.(k2)


#  Find e-vector:
#  Start with old method (not isotropic) for now.
#  Pick e1 if kk not parallel to ee1. ee2 else.
      if ky==0. && kz==0.
        ex = 0.; ey = 1.; ez = 0.
      else
        ex = 1.; ey = 0.; ez = 0.
      end

#  Isotropize ee in the plane perp. to kk by
#  (1) constructing two basis vectors for the plane perpendicular
#      to kk, and
#  (2) choosing a random direction in that plane (angle phi)
#  Need to do this in order for the forcing to be isotropic.

    kk = [kx, ky, kz]
    kk = collect(Iterators.flatten(kk))
    ee = [ex, ey, ez]
    e1 = cross(kk, ee)
    norm_e1 = norm(e1)
    e1 = e1 / norm(e1) 

    #call dot2(e1,norm); e1=e1/sqrt(norm) # e1: unit vector perp. to kk
#    call cross(kk,e1,e2)
#    call dot2(e2,norm); e2=e2/sqrt(norm) # e2: unit vector perp. to kk, e1
    e2 = cross(kk, e1)
    e2 = e2 / norm(e2) 
   # call random_number_wrapper(phi,CHANNEL=channel_force); phi = phi*2*pi
    phi = pi*2*rand(1)

    ee = cos.(phi) .* e1 .+ sin.(phi) .* e2
    ex = ee[1]; ey = ee[2]; ez = ee[3]
    #
#  k.e
#
      kde = dot(kk,ee)
#
#  k x e
#
      kex=ky*ez-kz*ey
      key=kz*ex-kx*ez
      kez=kx*ey-ky*ex
#
#  k x (k x e)
#
      kkex = ky .* kez - kz .* key
      kkey = kz .* kex - kx .* kez
      kkez = kx .* key - ky .* kex
#
##  ik x (k x e) + i*phase
#
#  Normalize ff; since we don't know dt yet, we finalize this
#  within timestep where dt is determined and broadcast.
#
#  This does already include the new sqrt(2) factor (missing in B01).
#  So, in order to reproduce the 0.1 factor mentioned in B01
#  we have to set force=0.07.
#
#  Furthermore, for |relhel| < 1, sqrt(2) should be replaced by
#  sqrt(1.+relhel**2). This is done now (9-nov-02).
#  This means that the previous value of force=0.07 (for relhel=0)
#  should now be replaced by 0.05.
#
#  Note: kav is not to be scaled with k1_ff (forcing should remain
#  unaffected when changing k1_ff).
#
      relhel = 0. # "sigma" in Eq.227 of the menu 
      cs0eff = 1.
      slope_ff = 0.
#      ffnorm = sqrt(1. + relhel^2) *k*sqrt(k2-kde^2)/sqrt(kav*cs0eff^3)*(k/kav)^slope_ff
ffnorm = sqrt.(1. .+ relhel .^ 2) .* k .* sqrt.(abs.(k2 .- kde .^ 2)) ./ sqrt.(kav .* cs0eff.^ 3) .* (k ./ kav) .^ slope_ff

#
#  need to multiply by dt (for Euler step), but it also needs to be
#  divided by sqrt(dt), because square of forcing is proportional
#  to a delta function of the time difference
#

      dt_force = dt * dt_force_over_lbm
      fact = force_factor ./ ffnorm .* sqrt.(dt_force) # need to specify "dt"

#
#  prefactor; treat real and imaginary parts separately (coef1 and coef2),
#  so they can be multiplied by different profiles below.
#
    coef1 = k .* [kex,key,kez]
    coef2= relhel * [kkex,kkey,kkez]
#
#   return  coef1, coef2, kk, phase, fact, fda
   return  coef1, coef2, kk, phase, fact
end

#function fxyz_coefs_hel(x, y, z, coef1,coef2,coef3,kk,phase,fact) 
function fxyz_coefs_hel(x, y, z, kk, phase, fact) 
    k1_ff = 1.
fx = exp.((kk[1] .* k1_ff .* x .+ phase)*im) .* fact # N*exp{} of Eq.226 in the manual. Note that the "N" and "phi" terms are included in the "fx" 
    fy = exp.((kk[2] .* k1_ff .* y)*im)
    fz = exp.((kk[3] .* k1_ff .* z)*im)
    #
    return fx,fy,fz
end    

## Forcing for Poiseuille plane flow:
#function force_Poiseuille(body_force_Poiseuille)
#    fx_domain = zeros(LX, LY) 
#    fx_domain .= body_force_Poiseuille
#    fy_domain = zeros(LX, LY)
#
#    return fx_domain, fy_domain
#end
## Forcing for Poiseuille plane flow.
