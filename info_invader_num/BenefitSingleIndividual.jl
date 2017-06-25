# expected fitness of an informed mutant (0 < I <=1) invading a resident
# density-independently emigrating population

# IMPORTANT NOTICE: In modern versions of Julia it is sometimes necessary to use
# a dot for mathematical operations that shall be applied to all elements of an array;
# e.g. if x = [1,2,3,4,5] and y = 5: z = x .* 5 instead of x * 5

using DataFrames

popdata = readtable("base_model_noinvasion/popsizes.wsv")

Start = time()
NEVAL = 1000 # the actual number of used population sizes
preN = popdata[1:NEVAL,:pre_disp] # we extract only these
postN = popdata[1:NEVAL,:post_disp]
Ntot = sum(preN) # total population size
# simulation parameters
K = 100
μ = 0.1
Λ = 2.0
σ = 1.0

# degrees of information to be tested
I = linspace(0,1,21)

# emigration thresholds to be tested (77 means steps of size 5)
vCT = linspace(20,400,77)

α = (Λ - 1.) / K # term α in Beverton-Holt equation
exp_fert = Λ ./ (1.0 .+ postN .* α) # expected payoff in each patch AFTER dispersal
mfit_emi = mean(exp_fert) # expected mean payoff of successful dispersers
mfit_emimu = mfit_emi * (1. - μ) # the same, but including dispersal

emi = zeros(Float64,length(vCT),length(I))
gain = zeros(Float64,length(vCT),length(I))
#sum(knowledge[1]*rancomp.>vCT[k])/length(knowledge)

for i in 1:length(I)
  println("I=",I)
  rancomp = (1.-I[i]) .* preN
  knowledge = I[i] .* preN
  for k in 1:length(vCT)
    println("  CV = ",vCT[k])
    flush(STDOUT)
    if I==1
      pemiNBCT = convert(Int64,preN.>vCT[k])
    else
      pemiNBCT = map(x->sum(x+rancomp.>vCT[k]),knowledge)./length(knowledge) # the mapping means basically applying the function x to each element of knowledge
    end
    
    fitnessNBCT = pemiNBCT .* mfit_emimu .+ (1.0.-pemiNBCT) .* exp_fert
    emi[k,i] = sum(preN .* pemiNBCT)/Ntot
    gain[k,i] = sum(preN .* fitnessNBCT)/Ntot
  end
end

End = time()

print(End-Start)

writedlm("gain.wsv",gain,'\t')
writedlm("emi.wsv",emi,'\t')
