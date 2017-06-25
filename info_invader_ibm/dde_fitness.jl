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

const REPLI_START = parse(Int64,ARGS[1]) # performed replicates from command line
const REPLI_END = parse(Int64,ARGS[2])
const II = parse(Float64,ARGS[3])

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
  d::Float64 # trait of density-independent dispersal
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
  d::Array{Float64,1}
  d_sd::Array{Float64,1}
  cth_t::Array{Float64,1}
  cth_t_sd::Array{Float64,1}
  meta_n::Array{Int64,1}
  expressed_Cth::Array{Float64,1}
  ω::Array{Float64,1}
  invader_b4::Int64 # number of dde individuals before dispersal
  invader_disp::Int64 # number of dde individuals dispersed
  residents_b4::Int64 # same for die individuals
  residents_disp::Int64
  memory::Array{Float64,2}
  N_guess::Float64
  adults_t::Int64
  kids_t::Int64
  emi_invader::Array{Float64,1}
  emi_residents::Array{Float64,1}
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
        ind = TInd(0.5+rand(),rand(),false,true)
      else
        # create a TAe individual
        ind = TInd(0.5+rand(),rand(),false,false)
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
  d = 0.
  if ind.T
    if (ana.N_guess/par.K)>=ind.Cth_T
      d = 1.
    end
  else
      d = ind.d
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
  ana.invader_b4 = 0
  ana.invader_disp = 0
  ana.residents_b4 = 0
  ana.residents_disp = 0
  ana.adults_t = 0

  # for each patch
  for p = 1:NPATCH

    i = 0

    while i<size(world[p].ind,1)

      i += 1

      if !world[p].ind[i].disp

        ana.N_guess = guess_N(p,ana,par)

        # count individuals before dispersal
        world[p].ind[i].T ? ana.invader_b4 += 1 : ana.residents_b4 += 1
        ana.adults_t = world[p].ind[i].T ? ana.adults_t + 1 : ana.adults_t

        P_d = disp_prob(world[p].ind[i],ana,p,par)

        # if dispersal
        if rand()<P_d

          # count emigrants
          world[p].ind[i].T ? ana.invader_disp += 1 : ana.residents_disp += 1

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

  newborn = TInd(parent.Cth_T,parent.d,false,parent.T)

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
  ds = Float64[]

  for p = 1:NPATCH # for every patch
    pops = size(world[p].ind,1) # get popsize
    metan = metan + pops
    for i = 1:pops
      push!(ds,world[p].ind[i].d) # store alleles
      push!(cths_t,world[p].ind[i].Cth_T)
    end
  end

  push!(ana.t,t)
  push!(ana.meta_n,metan)
  if metan>0
    push!(ana.d,mean(ds))
    push!(ana.d_sd,std(ds));
    push!(ana.cth_t,mean(cths_t))
    push!(ana.cth_t_sd,std(cths_t))
    push!(ana.emi_residents,ana.residents_disp/ana.residents_b4)
    push!(ana.emi_invader,ana.invader_disp/ana.invader_b4)
  else
    push!(ana.d,0.0)
    push!(ana.d_sd,0.0)
    push!(ana.cth_t,0.0)
    push!(ana.cth_t_sd,0.0)
    push!(ana.emi_residents,0.0)
    push!(ana.emi_invader,0.0)
  end

  println(  outfile,t,"  ",metan,"  ",ana.residents_disp/ana.residents_b4,"  ",
            ana.invader_disp/ana.invader_b4,"  ",ana.d[end],"  ",ana.d_sd[end],
            "  ",ana.cth_t[end],"  ",ana.cth_t_sd[end],"  ",ana.ω[end])
  flush(outfile) # to make sure that the buffer is directly transferred

end

function plot(  t::Int64,
                ana::TAnalysis  )
  emi_inv = ana.emi_invader
  emi_res = ana.emi_residents
  cth = ana.cth_t
  ts = ana.t
  gain = ana.ω
  @rput ts emi_inv cth emi_res gain ANAMOD

  R"gain = gain[seq(ANAMOD,length(gain),ANAMOD)]"
  R"plot(emi_res~ts,type='b',lwd=2,pch=16,bty='l',xlab='time',ylab='emi',ylim=range(c(emi_res,emi_inv)))"
  R"lines(emi_inv~ts,type='b',lwd=2,pch=1,col='grey75')"
  R"legend('topright',lwd=2,pch=c(16,1),col=c('black','grey75'),c('resident','invader'),cex=1.75)"
  R"plot(cth~ts,type='b',lwd=2,pch=16,bty='l',xlab='time',ylab='threshold density (invader)')"
  R"plot(gain~ts,type='b',lwd=2,pch=16,bty='l',xlab='time',ylab='fitness gain')"
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
  ana = TAnalysis(  Int64[],Float64[],Float64[],Float64[],Float64[],
                    Int64[],Float64[],Float64[],0,0,0,0,
                    zeros(MEMORY,NPATCH),0.,0,0,Float64[],Float64[] )

  # read in parameters
  par = TParameters(  f[1,:K],f[1,:μ],f[1,:λ],
                      f[1,:σ],f[1,:tmax],II  )

  # perform replicate simulations
  for r = REPLI_START:REPLI_END

    # write popsizes only in the first run
    WRITE_POPSIZES = r==1 ? true : false

    # initialize output file
    outfilename = join(["results/output_I",par.I,"_rep",r,".txt"])
    outfile = open(outfilename,"w")
    println(outfile,"t  N  emi_res  emi_inv  d  d_sd  Cth_T  Cth_T_sd  gain")

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
      println("t = ",t)

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
      if t%ANAMOD == 0 #&& t > par.tmax
        analyze(t,outfile,world,ana)
 #       plot(t,ana)
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

#R"x11(width=14,height=5)"
#R"par(mfrow=c(1,3),cex.lab=2,cex.axis=1.75,mar=c(5,5,1,1))"
run_package()
