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
using RCall

#######################
#### const declarations

# the number of chosen individuals to follow the alternative strategy
# relative to total N:
const THRESH_T = 10
# the number of generations over which population densities are stored for
# randomly drawn population sizes:
const MEMORY = 5
# number of patches:
NPATCH = 1000
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
  d::Float64 # trait of density-independent dispersal
  disp::Bool  # has the individual dispersed already?
end

# a patch is just defined by the individuals inhabiting it
type TPatch
  ind::Array{TInd,1}
end

# necessary simulation parameters
type TParameters
  K::Float64    # carrying capacity
  μ::Float64    # dispersal mortality
  λ::Float64    # per capita growth rate
  σ::Float64    # degree of environmental stochasticity
  tmax::Int64   # no. of simulated generations
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
      ind = TInd(rand(),false)
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


# the actual dispersal function: ask each individual, whether it will emigrate,
# delete it from its natal patch and, if it survives transition, put it into
# a randomly determined target patch
function dispersal( world::Array{TPatch,1},
                    par::TParameters  )


  # for each patch
  for p = 1:NPATCH

    i = 0

    while i<size(world[p].ind,1)

      i += 1

      if !world[p].ind[i].disp


        # if dispersal
        if rand()<world[p].ind[i].d

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

  newborn = TInd(parent.d,false)

  # now with decreasing mutation rates for the TAe strategy
  mut = MUTSTART * exp( -MUTDEC * t/par.tmax)

  if (rand()<mut)
    newborn.d += 0.2*(rand()-0.5)
  end

  # making clear that d lies in the interval [0,1]
  newborn.d = newborn.d < 0. ? 0. :
              newborn.d > 1. ? 1. :
              newborn.d

  return newborn
end


# the actual reproduction function: for each patch density-dependent survival
# probability s is calculated, based on Beverton & Holt (1957), each individual
# gives birth to a number of offspring drawn from a Poisson distribution and
# it is stored, whether the kids have been given birth by a T emigrant or not;
# in the end these values are used to calculate relative fitness ω of T
function reproduction(  t::Int64,
                        world::Array{TPatch,1},
                        par::TParameters  )

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

end

##############################
### analysis-related functions

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


##################################
### simulation-related function(s)

# the actual simulation is performed by this function
function simulate(  ex::Int64,
                    f::DataFrames.DataFrame  )

  # initialize variables
  world = TPatch[]

  # read in parameters
  par = TParameters(  f[1,:K],f[1,:μ],f[1,:λ],
                      f[1,:σ],f[1,:tmax]  )

  # perform replicate simulations

    # initialize array for storing population sizes before and after disp.
    # and the according file with appropriate filenames
    popsizes = zeros(Int64,NPATCH,2)
    popoutfile = open("popsizes.wsv","w")
    println(popoutfile,"t  pre_disp  post_disp")


    # initialize world
    init_world(world,par)

    # for each generation (show the progress in the terminal, important is only the loop)
    for t = 1:(par.tmax)
      println("t = ",t)

      # 1. let individuals disperse (and probably store & write population sizes)
      if (t>9000)
        store_popsizes(t,popsizes,world,true,par)
      end

      dispersal(world,par)

      if (t>9000)
        store_popsizes(t,popsizes,world,false,par)
        write_popsizes(t,popsizes,popoutfile,par)
      end

      # 2. let individuals reproduce
      reproduction(t,world,par)


    end

    close(popoutfile)
end


# for each line in parameters.csv, a simulation package with optional
# replicates is started
function run_package()

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
