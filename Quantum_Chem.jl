using Fermi
using Fermi.Integrals
using Plots
using ProgressMeter
using Molecules
using Unitful
using PhysicalConstants.CODATA2018:a_0,e,ε_0,ħ,c_0,α,h,k_B

## This section calculates the potential energy curve for the ground state
# of the CaH system (as in the Habli et al Paper)

# Distance vector for system
Rval2 = round.([0.99 + 0.01*i for i = 1:901],digits=2)
# results vector
E_rhf = []

@showprogress for r in Rval2
  # Set the CaH system
  mol = string("""
  Ca          0.0 0.0 0.0
  H           $r  0.0    0.0
  """);

  # settings for solve
  @set {
    basis cc-pvtz # coupled cluster doublet basis
    scf_max_iter	100 # iter max out
    DIIS true # convergence tool
    charge -1 # Charge of the system (free e⁻)
    multiplicity 3 # multiplicity of the molecule (2n+1)
  }

  Fermi.Options.set("molstring", mol);
  
  # Calc the energy at a given radius
  rhf_e = @energy rhf;

  # Save result, noting in atomic units
  push!(E_rhf, rhf_e.energy)

end

# Energy in atomic units to eV
E_h = (ħ*c_0*α)/(a_0) |> u"eV"

# Add units to the radius
new_r = (Rval2)u"angstrom"

# Scale change to indicate attraction/repulsion 
E_CaH = (E_rhf).*E_h .|> u"keV"#.-E_rhf[end]).*E_h

# plot
plot(new_r,E_CaH,label="Reduced Hartree Fock CaH X²Σ⁺",xlabel="R",ylabel="E")




## This section considers the Ca H₂ system which seems to be case in Proxima cen
rhf_min_CaH,rhf_idx = findmin(E_CaH)
rhf_sep = new_r[rhf_idx]


plot!(new_r[1:end-1],diff(E_CaH))
Rval = round.([0.9 + 0.01*i for i = 1:191],digits=2)
E_H_2 = []
@reset
@showprogress for r in Rval2
  # Set the CaH system
  mol = string("""
  Ca          0.0   0.0     0.0
  H           $r    0.36    0.0
  H           $r   -0.36    0.0
  """);

  # settings for solve
  @set {
    basis cc-pvtz
    scf_max_iter	100
    DIIS true # helper for convergence
    reference rhf #Restricted Hartree-Fock
    charge 0 # Charge of the system
    multiplicity 1 # multiplicity of the molecule
  }
  Fermi.Options.set("molstring", mol);
  #I = Fermi.Integrals.IntegralHelper(eri_type=Fermi.Integrals.Chonky())
  # Calc the energy at a given radius
  wf = @energy rhf;
  # Save result, noting E_h is in Hartrees
  push!(E_H_2, wf.energy)
end

#new_r = (Rval[1:end]u"angstrom")./a_0 .|> upreferred
E_CaH2 = (E_H_2[1:end]).*E_h .|> u"keV"#- E_H_2[end])*E_h
plot!(new_r,E_CaH2,label="CaH₂",xlabel="R(a.u.)")


new_test = (E_CaH .- E_CaH2) .|> u"eV"

spectrum = c_0./(new_test./h) .|> u"nm"

plot(spectrum)







@reset
@molecule {
  Ca 0 0 0
}
@set {
  basis cc-pvtz
  charge -2
  multiplicity 5
}

Ca_single = @energy rhf
Ca_single.energy










## Data points from the Habli paper
curves = [[17.454 -0.000068182],
[17.454 0.10077],
[17.454 0.10670]]

# distance between points out at 10 a.u.
dist = abs(curves[3][2]-curves[1][2])

# Energy in Hartrees
E_h = (ħ*c_0*α)/(a_0) |> u"eV"


dist = 0.0940738008715698
# Line core for CaH 
λ_habli = c_0*h/((dist*E_h)) |> u"nm"

c_0/((2.9u"eV")/h) |> u"nm"

# Ca I normal line core (NIST Transitional database)
λ=422.673u"nm"

# distance check using NIST value
dist_1 = c_0*h/(λ*E_h) |> upreferred

# This provides the line core, lorentzian wings from the sum of the transition probabilities


## Taking the EOS code that Mike Ireland developed to understand Debye shielding and interatomic distances
T = ([1500.0, 1660.0, 1930.0, 2350, 3320])u"K"

abund_H = [-0.24137298  0.8650198   2.36268973  3.96646332  5.95206568] |> vec
abund_Ca =  [-2.31584277 -1.61687306 -0.61687873  0.38295317  1.35964597] |> vec
abund_H2 = [3.04282275 3.7416172  4.74099784 5.7382488  6.70517299] |> vec

# Partial Pressures of HI, Ca and H2 in Proxima Cen like atmosphere
p_HI = 10 .^(abund_H).*1u"dyn/cm^2"
p_Ca = 10 .^(abund_Ca).*1u"dyn/cm^2"
p_H2 = 10 .^(abund_H2).*1u"dyn/cm^2"

# Conversion to number density
n_HI = p_HI./(k_B*T) .|> u"cm^(-3)"
n_Ca = p_Ca./(k_B*T) .|> u"cm^(-3)"
n_H2 = p_H2./(k_B*T) .|> u"cm^(-3)"

# This is the debye length
D_HI = sqrt.(ε_0*k_B .*T ./(n_HI .*e^2)) .|> u"angstrom"
D_H2 = sqrt.(ε_0*k_B .*T ./(n_H2 .*e^2)) .|> u"angstrom"
D_Ca = sqrt.(ε_0*k_B .*T ./(n_Ca .*e^2)) .|> u"m"

# number of Ca particles in Debye Sphere
δ_HI = 4π/3 .*D_Ca.^3 .*n_HI .|> upreferred
δ_H2 = 4π/3 .*D_Ca.^3 .*n_H2 .|> upreferred
δ_Ca = 4π/3 .*D_Ca.^3 .*n_Ca .|> upreferred

# Mean interparticle distance
r_0_H2 = (4π .*n_H2./3).^(-1/3) .|> u"angstrom"
r_0_HI = (4π .*n_HI./3).^(-1/3) .|> u"angstrom"
r_0_Ca = (4π .*n_Ca./3).^(-1/3) .|> u"angstrom"

# number of particles in mean Ca interparticle sphere
n_s_HI = 4π/3 .*r_0_Ca.^3 .*n_HI .|> upreferred
n_s_H2 = 4π/3 .*r_0_Ca.^3 .*n_H2 .|> upreferred
n_s_Ca = 4π/3 .*r_0_Ca.^3 .*n_Ca .|> upreferred


# This means that we expect to see Calcium particles within out Debye sphere
@assert r_0_H2 < D_H2

# Range of interaction


β = (r_0_HI[1]./r).^2 .|> upreferred

W(x) = 3/2*x^(-5/2)*exp(-x^(-3/2))

Linear_Stark = W.(β)


# Probability function of H2 in range (r,r+dr), not including Debye shielding
W(r,N) = 4π*r^2*N*exp(-4/3*π*r^3*N)

using Plots
r = (0.01:0.1:60)u"angstrom"
sphere_prob = W.(r,n_H2[1]) .|> u"cm^(-1)"
sphere_prob1 = W.(r,n_H2[2]) .|> u"cm^(-1)"
sphere_prob2 = W.(r,n_H2[3]) .|> u"cm^(-1)"
sphere_prob3 = W.(r,n_H2[4]) .|> u"cm^(-1)"
sphere_prob4 = W.(r,n_H2[5]) .|> u"cm^(-1)"
sphere_prob = sphere_prob./sum(sphere_prob)
sphere_prob1 = sphere_prob1./sum(sphere_prob1)
sphere_prob2 = sphere_prob2./sum(sphere_prob2)
sphere_prob3 = sphere_prob3./sum(sphere_prob3)
sphere_prob4 = sphere_prob4./sum(sphere_prob4)
plot(r,sphere_prob)
plot!(r,sphere_prob1)
plot!(r,sphere_prob2)
plot!(r,sphere_prob3)
plot(r,sphere_prob4)




# now considering the Holtzmark approach
using QuadGK
p=6
W_H(β) = 2*β/π*quadgk(y -> exp(-y^(3/p))*y*sin(β*y),0,Inf)[1]


Holtz = W_H.(β)

n=3
plot(r[n:end],Holtz[n:end],label="Holtzmark Distribution based of Debye Shielding")
plot!(r[n:end],Linear_Stark[n:end],label="Linear Stark Effect")