#############################################################
### The fitness value of density-dependent dispersal
### - source code of the simulations -
#############################################################

# Dr. Alexander Kubisch
# Theoretical Evolutionary Ecology Group
# Department of Animal Ecology and Tropical Biology
# University of Würzburg
# Germany
# alexander.kubisch@uni-wuerzburg.de
# https://github.com/akubisch
# 2016
# MIT License

##################################
### using some additional packages

using Distributions # provides statistical distributions (incl. random generation)
using DataFrames # for reading in data

#######################
#### const declarations

const REPLI_START = parse(Int64,ARGS[1]) # performed replicates from command line
const REPLI_END = parse(Int64,ARGS[2])

# the number of chosen individuals to follow the alternative strategy
# relative to total N:
const N_CHOSEN = 0.01
const THRESH_T = 10
# the number of generations over which population densities are stored for
# randomly drawn population sizes:
const MEMORY = 5
# number of patches:
NPATCH = 1000
# interval for analysis:
ANAMOD = 10
# additional simulation time for analysis
ANATIME = 1000
# initial mutation rate
MUTSTART = 0.001
# decrease rate of mutation rate relative to t/par.tmax
MUTDEC = 0.0

####################
### type definitions

# an individual is defined by its density-dependent dispersal threshold
# according to the extended asymptotic strategy (Cth_TAe),
# its "simple" dispersal threshold Cth_T, a marker, whether it has
# already dispersed and another marker defining, whether it follows
# the simple threshold model (i.e. is an invader) or the TAe (resident)
type TInd
  Cth_T::Float64 # trait of simple threshold model
  Cth_TAe::Float64 # trait of asymptotic extended threshold model
  γ::Float64    # parameter shaping the extended threshold
  disp::Bool  # has the individual dispersed already?
  T::Bool # is this individual fully informed with perfect strategy?
end

# a patch is just defined by the individuals inhabiting it
type TPatch
  ind::Array{TInd,1}
end

# some simple analyses of our results
type TAnalysis
  t::Array{Int64,1}
  cth_tae::Array{Float64,1}
  cth_tae_sd::Array{Float64,1}
  cth_t::Array{Float64,1}
  cth_t_sd::Array{Float64,1}
  γ::Array{Float64,1}
  γ_sd::Array{Float64,1}
  meta_n::Array{Int64,1}
  expressed_Cth::Array{Float64,1}
  ω::Array{Float64,1}
  t_b4::Int64 # number of dde individuals before dispersal
  t_disp::Int64 # number of dde individuals dispersed
  tae_b4::Int64 # same for die individuals
  tae_disp::Int64
  memory::Array{Float64,2}
  N_guess::Float64
  adults_t::Int64
  kids_t::Int64
end

# necessary simulation parameters
type TParameters
  K::Float64    # carrying capacity
  μ::Float64    # dispersal mortality
  λ::Float64    # per capita growth rate
  σ::Float64    # degree of environmental stochasticity
  tmax::Int64   # no. of simulated generations
  I::Float64    # accuracy of resident information
end

######################################
### basic functions and initialization

# log-normal random number generation based on the mean and s.d. of the
# resulting distribution
function logn(  μ::Float64,
                s_d::Float64  )
  σ = sqrt(log((s_d/μ)^2 + 1))
  ζ = log(μ^2/sqrt(s_d^2 + μ^2))
  distri = LogNormal(ζ,σ)
  return rand(distri)
end

# initialize the world: read in parameters, populate all patches with indiv.
# assign random traits to these and define a fraction of N_CHOSEN to perform
# density-dependent emigration
function init_world(  world::Array{TPatch,1},
                      par::TParameters  )

  for p = 1:NPATCH # for every patch

    inds = TInd[] # empty container for individuals

    for i = 1:par.K # K times
      if rand() < N_CHOSEN
        # create a T individual
        ind = TInd(10.0*rand(),0.5+rand(),rand(),false,true)
      else
        # create a TAe individual
        ind = TInd(10.0*rand(),0.5+rand(),rand(),false,false)
      end
      # and put it into the container
      push!(inds,ind)
    end

    patch = TPatch(inds) # transform the container into the TPatch form
    push!(world,patch) # and put it into our world
  end

end

###############################
### dispersal-related functions

# define a target patch for dispersal - must be different from natal
function find_target_patch( p::Int64 )

  target = rand(1:NPATCH)

  while target==p
    target = rand(1:NPATCH)
  end

  return target

end

# this function returns the emigration probability based on the individual's
# traits and, respective to its strategy, population density before dispersal
function disp_prob( ind::TInd,
                    ana::TAnalysis,
                    p::Int64,
                    par::TParameters  )
  d = 0
  if ind.T
    if ana.memory[MEMORY,p]/par.K>=ind.Cth_T
      d = 1
    end
  else
    if ana.N_guess/par.K>=ind.Cth_TAe
      d = 1. - (ind.Cth_TAe*(1.0-ind.γ))/(ana.N_guess/par.K-ind.Cth_TAe*ind.γ)
    end
  end
  return(d)
end

# the following function returns the estimated population density based
# on the degree of information accuracy (par.I)
function guess_N( p::Int64,
                  ana::TAnalysis,
                  par::TParameters  )
  if par.I > 0.99
    N_guess = ana.memory[MEMORY,p]
  else
    r = rand(0:(NPATCH*MEMORY-1))
    rpat = mod(r,NPATCH)+1
    rmem = div(r,NPATCH)+1
    N_rand = ana.memory[rmem,rpat]
    N_guess = par.I * ana.memory[MEMORY,p] + (1.-par.I)*N_rand
  end

  return N_guess
end

# the actual dispersal function: ask each individual, whether it will emigrate,
# delete it from its natal patch and, if it survives transition, put it into
# a randomly determined target patch
function dispersal( world::Array{TPatch,1},
                    ana::TAnalysis,
                    par::TParameters  )

  # resetting analysis
  ana.t_b4 = 0
  ana.t_disp = 0
  ana.tae_b4 = 0
  ana.tae_disp = 0
  ana.adults_t = 0

  # for each patch
  for p = 1:NPATCH

    i = 0

    while i<size(world[p].ind,1)

      i += 1

      if !world[p].ind[i].disp

        ana.N_guess = guess_N(p,ana,par)

        # count individuals before dispersal
        world[p].ind[i].T ? ana.t_b4 += 1 : ana.tae_b4 += 1
        ana.adults_t = world[p].ind[i].T ? ana.adults_t + 1 : ana.adults_t

        P_d = disp_prob(world[p].ind[i],ana,p,par)

        # if dispersal
        if rand()<P_d

          # count emigrants
          world[p].ind[i].T ? ana.t_disp += 1 : ana.tae_disp += 1

          # if survival
          if rand()>par.μ
            world[p].ind[i].disp = true
            # find a new patch
            newpatch = find_target_patch(p)
            # transfer the individual
            push!(world[newpatch].ind,world[p].ind[i])
          end

          # delete it from its natal patch
          deleteat!(world[p].ind,i)
          i -= 1
        end
      end
    end
  end

end

###################################
### reproduction-related procedures

# this function returns a newborn individual inheriting its parent's values;
# NOTE: while a newborn inherits Cth_T from its parent, this will later
# be redrawn
function genetics(  t::Int64,
                    parent::TInd,
                    par::TParameters  )

  newborn = TInd(parent.Cth_T,parent.Cth_TAe,parent.γ,false,parent.T)

  # now with decreasing mutation rates for the TAe strategy
  mut = MUTSTART * exp( -MUTDEC * t/par.tmax)

  if (rand()<mut)
    newborn.Cth_TAe += 0.2*(rand()-0.5)
  end
  if (rand()<mut)
    newborn.γ += 0.2*(rand()-0.5)
  end

  # making clear that γ lies in the interval [0,1]
  newborn.γ = newborn.γ < 0. ? 0. :
              newborn.γ > 1. ? 1. :
              newborn.γ

  newborn.Cth_TAe = newborn.Cth_TAe < 0. ? 0. :
                    newborn.Cth_TAe

  return newborn
end

# the following function collects and stores all genotypes of Cth that have
# been expressed (i.e. used) in the last generation; these will later be
# redistributed over the whole metapopulation to increase selection pressure
function collect_expressed_Cth( world::Array{TPatch,1} )

  # first we reset the expressed_Cth vector to be empty
  expr = Float64[]
  # now we fill it with all expressed Cth_T values
  for p = 1:NPATCH
    for i = 1:size(world[p].ind,1)
      if world[p].ind[i].T
        push!(expr,world[p].ind[i].Cth_T)
      end
    end
  end

  return expr
end

# this function actually redistributes the expressed values of Cth (and only
# these) across the metapopulation and mutates them with probability par.m
function redistribute_Cth(  t::Int64,
                            world::Array{TPatch,1},
                            ana::TAnalysis,
                            par::TParameters  )

  ana.expressed_Cth = collect_expressed_Cth(world)

  for p = 1:NPATCH
    for i = 1:size(world[p].ind,1)
      world[p].ind[i].Cth_T = rand(ana.expressed_Cth)
      mut = MUTSTART * exp( -MUTDEC * t/par.tmax)
      if (rand()<mut)
        world[p].ind[i].Cth_T += 0.2*(rand()-0.5)
        world[p].ind[i].Cth_T = world[p].ind[i].Cth_T < 0. ? 0. :
                                world[p].ind[i].Cth_T
      end
      # redefine the chosen ones
      world[p].ind[i].T = rand() < N_CHOSEN ? true : false
    end
  end
end

# the actual reproduction function: for each patch density-dependent survival
# probability s is calculated, based on Beverton & Holt (1957), each individual
# gives birth to a number of offspring drawn from a Poisson distribution and
# it is stored, whether the kids have been given birth by a T emigrant or not;
# in the end these values are used to calculate relative fitness ω of T
function reproduction(  t::Int64,
                        world::Array{TPatch,1},
                        ana::TAnalysis,
                        par::TParameters  )

  ana.kids_t = 0

  # for each patch
  for p = 1:NPATCH

    kids = TInd[]

    # get popsize
    pops = size(world[p].ind,1)

    # environmental stochasticity based on σ
    Λ = par.λ
    if par.σ>0
      Λ = logn(par.λ,par.σ)
    end

    # determine survival probability
    α = (par.λ - 1.0) / par.K
    s = 1.0 / (1.0 + α * pops)

    # creates the kid number distribution for this patch and year
    distri_kids = Poisson(Λ * s)

    # for each individual
    for i = 1:pops
      # draw number of kids from the dist'n
      no_kids = rand(distri_kids)
      # store the fecundity
      ana.kids_t = world[p].ind[i].T ? ana.kids_t + no_kids : ana.kids_t

      # for every kid
      for k = 1:no_kids
        newborn = genetics(t,world[p].ind[i],par)
        push!(kids,newborn)
      end

    end

    # adults get replaced by offspring
    world[p].ind = TInd[]
    for i = 1:size(kids,1)
      push!(world[p].ind,kids[i])
    end

  end

  # mean fitness is calculated and attached to the resp. analysis array
  push!(ana.ω,(ana.kids_t-ana.adults_t)/ana.adults_t)

end

##############################
### analysis-related functions

# compute mean trait values and population size and write these incl. fitness
# into a file
function analyze( t::Int64,
                  outfile::IOStream,
                  world::Array{TPatch,1},
                  ana::TAnalysis  )

  cths_t = Float64[] # container for all alleles in the metapop
  metan = 0 # metapopsize
  cths_tae = Float64[]
  γs = Float64[]

  for p = 1:NPATCH # for every patch
    pops = size(world[p].ind,1) # get popsize
    metan = metan + pops
    for i = 1:pops
      push!(cths_tae,world[p].ind[i].Cth_TAe) # store alleles
      push!(cths_t,world[p].ind[i].Cth_T)
      push!(γs,world[p].ind[i].γ)
    end
  end

  push!(ana.t,t)
  push!(ana.meta_n,metan)
  if metan>0
    push!(ana.cth_tae,mean(cths_tae))
    push!(ana.cth_tae_sd,std(cths_tae));
    push!(ana.cth_t,mean(cths_t))
    push!(ana.cth_t_sd,std(cths_t))
    push!(ana.γ,mean(γs))
    push!(ana.γ_sd,std(γs))
  else
    push!(ana.cth_tae,0)
    push!(ana.cth_tae_sd,0)
    push!(ana.cth,0)
    push!(ana.cth_sd,0)
  end

  println(  outfile,t,"  ",metan,"  ",ana.tae_disp/ana.tae_b4,"  ",
            ana.t_disp/ana.t_b4,"  ",ana.cth_tae[end],"  ",ana.cth_tae_sd[end],
            "  ",ana.cth_t[end],"  ",ana.cth_t_sd[end],"  ",ana.γ[end],"  ",
            ana.γ_sd[end],"  ",ana.ω[end])
  flush(outfile) # to make sure that the buffer is directly transferred
end

# the following function saves population sizes before or after dispersal
# if WRITE_POPSIZES is set to be true
function store_popsizes(  t::Int64,
                          pops::Array{Int64,2},
                          world::Array{TPatch,1},
                          pre::Bool,
                          par::TParameters  )
  # if it is late enough
  if t > par.tmax - THRESH_T
    # if this is before dispersal
    if pre
      for p in 1:NPATCH
        pops[p,1] = size(world[p].ind,1)
      end
    else
      for p in 1:NPATCH
        pops[p,2] = size(world[p].ind,1)
      end
    end
  end
end

# the following function saves these population sizes to a file
function write_popsizes(  t::Int64,
                          pops::Array{Int64,2},
                          outfile::IOStream,
                          par::TParameters  )
  if t > par.tmax - THRESH_T
    for p = 1:NPATCH
      println(outfile,t,"  ",pops[p,1],"  ",pops[p,2])
    end
  end
end

# the following function updates the individuals' "memory", including past
# and present population densities
function update_memory( world::Array{TPatch,1},
                        ana::TAnalysis  )
  for m = 1:(MEMORY-1)
    for p = 1:NPATCH
      ana.memory[m,p] = ana.memory[m+1,p]
    end
  end

  for p = 1:NPATCH
    ana.memory[MEMORY,p] = size(world[p].ind,1)
  end
end

##################################
### simulation-related function(s)

# the actual simulation is performed by this function
function simulate(  ex::Int64,
                    f::DataFrames.DataFrame  )

  # initialize variables
  world = TPatch[]
  ana = TAnalysis(  Int64[],Float64[],Float64[],Float64[],Float64[],Float64[],
                    Float64[],Int64[],Float64[],Float64[],0,0,0,0,
                    zeros(MEMORY,NPATCH),0.,0,0 )

  # read in parameters
  par = TParameters(  f[ex,:K],f[ex,:μ],f[ex,:λ],
                      f[ex,:σ],f[ex,:tmax],f[ex,:I]  )

  # perform replicate simulations
  for r = REPLI_START:REPLI_END

    # write popsizes only in the first run
    WRITE_POPSIZES = r==1 ? true : false

    # initialize output file
    outfilename = join(["results/output_exp",ex,"_rep",r,".txt"])
    outfile = open(outfilename,"w")
    println(outfile,"t  N  emi_TAe  emi_T  Cth_TAe  Cth_TAe_sd  Cth_T  Cth_T_sd  gamma  gamma_sd  gain")

    # initialize array for storing population sizes before and after disp.
    # and the according file with appropriate filenames
    if WRITE_POPSIZES
      popsizes = zeros(Int64,NPATCH,2)
      filename = join(["results/popsizes_I",par.I,".txt"])
      popoutfile = open(filename,"w")
      println(popoutfile,"t  pre_disp  post_disp")
    end

    # initialize world
    init_world(world,par)

    # for each generation (show the progress in the terminal, important is only the loop)
    for t = 1:(par.tmax+ANATIME)

      # 1. update stored population densities
      update_memory(world,ana)

      # 2. let individuals disperse (and probably store & write population sizes)
      if WRITE_POPSIZES
        store_popsizes(t,popsizes,world,true,par)
      end

      dispersal(world,ana,par)

      if WRITE_POPSIZES
        store_popsizes(t,popsizes,world,false,par)
        write_popsizes(t,popsizes,popoutfile,par)
      end

      # 3. let individuals reproduce
      reproduction(t,world,ana,par)

      # 4. collect and redistribute expressed genotypes (and define new chosen ones)
      redistribute_Cth(t,world,ana,par)

      # 5. analyze and save the results
      if t%ANAMOD == 0 && t > par.tmax
        analyze(t,outfile,world,ana)
      end

    end

    close(outfile)
    if WRITE_POPSIZES
      close(popoutfile)
    end
  end

end

# for each line in parameters.csv, a simulation package with optional
# replicates is started
function run_package()

  # create the results directory, if it is not existent already
  if !isdir("results")
    mkdir("results")
  end

  # read in parameters for each simulation run
  f = readtable("parameters.csv")

  # finally simulate each experiment
  for ex = 1:size(f,1)
    simulate(ex,f)
  end

end

# Last but not least: start the simulation!

println("----------------------------------------------------")
println("- The fitness value of density-dependent dispersal -")
println("----------------------------------------------------")
println()
println("by Alexander Kubisch, University of Würzburg")
println("alexander.kubisch@uni-wuerzburg.de")
println("https://github.com/akubisch/info_fitness")
println("MIT License")
println()
flush(STDOUT)

run_package()
