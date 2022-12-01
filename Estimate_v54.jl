# Constrained water and energy balance estimation 
# v4.2: adapted from v4.1 by treating SCF known, G estimated offline
# v4.3: adapted from v4.2 & 3.3 to handle to batch processing on OSC
# v4.4: adapted from v4.3 to update ipopt deprecation
# v4.5: moving σWRFG based on SCF 7/2/20 JLD
# v4.6: adding SDC to objective function: MD & JD: 9/28/20
# v4.7: adding "pseudo-valid prior" to starting poitns. remove SDC. MD: 1/7/21
# v4.8: tweak error parameters. MD & JD: 1/8/21
# v4.9 prior valid added to objective function. daily G capped in constraints. more outputs written 10/5/21
# v4.9.1 remove SWE from objective
# v4.9.2 added air temp
# v51 Converted to function but still uses textfile as input
# v52 Prototyping to pass the array directly to function [and fixed tab to 4 spaces]
# v53 To set tab to 4 spaces (by copying/pasting)
# v54 combine 9 separate textfile output in one txt file (due to file count limitation on Discover)


using JuMP
using Ipopt
using DelimitedFiles
# using Dates, CSV, DataFrames
# function blender(DataDir, exp_dir, WRFSWE, WRFP, WRFG, MSCF, AirT)
function blender(exp_dir, WRFSWE, WRFP, WRFG, MSCF, AirT)
    """
    Note: keyword argument defined after positional with ; if no default value provided, it is required
    Inputs
    ====== 
    exp_dir: Full path to exporting text outputs; created inside this script if it does not exist
    WRFSWE, WRFP, WRFG, MSCF, AirT : input variables

    Return : Nothing; only output text files created

    """
    # or pass v is a vector of arrays/matrix
    # exp_dir = "$DataDir/outputs2"  # BY
    # println("Ouput Dir:", exp_dir) # BY

    # # 0. read in data 
    # WRFSWE=readdlm(DataDir * "/WRFSWE.txt"); #[m]
    # WRFP=readdlm(DataDir * "/WRFP.txt"); #[m]
    # WRFG=readdlm(DataDir * "/WRFG.txt");
    # MSCF=readdlm(DataDir * "/MODSCAG.txt");
    # AirT=readdlm(DataDir * "/WRFT.txt");
    AirT=AirT.-273.15

    # # 1. define parameters
    # # 1.1 experiment parameters
    # exp_dir = DataDir 

    #σWSWE=0.2 #used for first run of 4.7
    σWSWE=0.4 # m
    RelPUnc=0.3 #[-]
    #σWRFG=25 # W m^-2 ... used for first run of 4.7
    σWRFG=15 # W m^-2
    σWMPmin=.001 #minimum uncertainty in fluxes: 1 mm/day
    σSCF=.1 # combined uncertainty in observed SCF and SDC
    nt = length(WRFSWE)  # 366  #365  # TODO this need to be updated based on lenght of input vector rather than fixed 365 days
    t=1:nt
    SWEmin=1.0e-6 #1/1000 mm
    ρnew=100 #density of new snow
    Δt=86400 #daily
    Gmax=300 # to prevent craziness. more than adequate for daily
    Gmin=-300 # to prevent craziness. more than adequate for daily
    ρnew=100 #density of new snow
    z0=0.01 # roughness length used in SDC
    mf=1.0 # melt factor used in SDC
    ρWRF=350 #note: previous versions used time-varying snow density but that is not currently exported

    # 1.1.1. read in σWRFG based on SCF

    #currently unused until further testing

    #σWRFG=zeros(nt,1)
    #for i=1:nt
    #  if MSCF[i]==0
    #	if WRFSWE[i]>0
    #	   σWRFG[i]=25;
    #	else 
    #	   σWRFG[i]=15;
    #       end
    #  else 
    #      σWRFG[i]=15;
    #  end
    #end

    # 1.2 physical parameters
    ρw=1000 #density of water
    Lf=0.334E6 #Latent heat of fusion J/kg
    cs_ice=2102 #specific heat capacity for ice J/kg/K
    
    # 1.3 Match up SWE and MSCF
    for i=1:nt
        if MSCF[i]==0  # Nov 04, 2022: "ERROR: LoadError: TypeError: non-boolean (Missing) used in boolean context" because some days there was no MODIS data
            WRFSWE[i] = 0
        end
    end
    
    # 2. compute useful variables and set up  arrays
    # 2.1 define SWEmax: upper limit of SWE, based on observed SCF
    SWEmax=zeros(nt,1)
    for i=1:nt
        if MSCF[i]==0
            SWEmax[i]=SWEmin
        else
            SWEmax[i]=5.0
        end
    end
    
    # 2.2 SDC denominator: dependent only on snow density
    DenomSDC=2.5.*z0.*(ρWRF./ρnew).^mf
    
    # 2.3 Uncertainty for accumulation 
    σWRFP=zeros(nt,1)
    for i=1:nt
        if WRFP[i]==0.
            σWRFP[i]=σWMPmin
        else
            σWRFP[i]=WRFP[i]*RelPUnc
        end
    end
    
    # 3. Solve
    # 3.1 Solve for G    
    Gmelt_pv=zeros(nt,1)
    G_pv=zeros(nt,1)
    U_pv=zeros(nt,1)
    SWEpv=WRFSWE
    for i=2:nt
        Gmelt_pv[i-1]=-(SWEpv[i]-SWEpv[i-1]-WRFP[i-1])*Lf*ρw/Δt
        if Gmelt_pv[i-1]<0.
            Gmelt_pv[i-1]=0
            SWEpv[i]=SWEpv[i-1]+WRFP[i-1]
        end
    end
    for i=2:nt
        if Gmelt_pv[i]>0. && MSCF[i] >0.
            G_pv[i]=Gmelt_pv[i]/MSCF[i]
            U_pv[i]=0.
        else
            G_pv[i]=WRFG[i]
            U_pv[i]=U_pv[i-1]+WRFP[i-1]*AirT[i-1]*ρw*cs_ice + G_pv[i-1]*MSCF[i-1]*Δt
            if U_pv[i]>0. || SWEpv[i]==0.
                U_pv[i]=0.
            end
        end
    end

    # 3.2 Solve for the posterior using prior valid
    m = Model(optimizer_with_attributes(Ipopt.Optimizer))
    @variable(m, SWEmin <= SWE[i=1:nt] <= SWEmax[i],start=SWEpv[i] )
    @variable(m,  Precip[i=1:nt]>=0. ,start=WRFP[i])
    @variable(m,  G[i=1:nt]<=Gmax , start=G_pv[i])
    @variable(m, Gmelt[i=1:nt]>=0, start=Gmelt_pv[i])
    @variable(m, Us[i=1:nt] <=0, start=U_pv[i])
    #@objective(m,Min,sum( (SWE-SWEpv).^2 ./σWSWE.^2)+sum((Precip-WRFP).^2 ./σWRFP.^2)+
    #                   sum((G-G_pv).^2 ./σWRFG.^2) )
    @objective(m,Min,sum((Precip-WRFP).^2 ./σWRFP.^2)+ sum((G-G_pv).^2 ./σWRFG.^2) )
    for i in 1:nt-1
        @constraint(m,SWE[i+1]==SWE[i]+Precip[i]-Gmelt[i]*Δt/Lf/ρw)
        @NLconstraint(m,Us[i+1]==Us[i]+(1-(tanh(Us[i]/10000)+1))*G[i]*MSCF[i]*Δt+
                                Precip[i]*AirT[i]*ρw*cs_ice) #m x K x kg/m3 x J/kg/K
    end
    @constraint(m,Us[1]==0) 
    @constraint(m,Us[nt]==0)
    for i in 1:nt
        @NLconstraint(m,Gmelt[i]==G[i]*MSCF[i]*(tanh(Us[i]/10000)+1))
    end
    
    # optimize!(m)
    # BY Create output folder; no error if the folder already exist
    mkpath(exp_dir)  # mkdir
    log_file =  "$exp_dir/log.txt"
    redirect_stdio(stdout=log_file, stderr=log_file) do
        optimize!(m)
        # println("Finish Optimation $DataDir")  
    end
    # 4. output
    SWEhat=JuMP.value.(SWE)
    GmeltHat=JuMP.value.(Gmelt)
    Ghat=JuMP.value.(G)
    Ushat=JuMP.value.(Us)
    Phat=JuMP.value.(Precip)

    out_vars = hcat(SWEhat, GmeltHat, Ghat, Phat, Ushat, G_pv, Gmelt_pv, U_pv, SWEpv)
    writedlm(exp_dir * "/out_vars.txt", out_vars)

    # writedlm(exp_dir * "/SWE.txt",SWEhat)
    # writedlm(exp_dir * "/Gmelt.txt",GmeltHat)
    # writedlm(exp_dir * "/G.txt",Ghat)
    # writedlm(exp_dir * "/Precip.txt",Phat)
    # writedlm(exp_dir * "/Us.txt",Ushat)
    # writedlm(exp_dir * "/Gpv.txt",G_pv)
    # writedlm(exp_dir * "/Gmeltpv.txt",Gmelt_pv)
    # writedlm(exp_dir * "/Upv.txt",U_pv)
    # writedlm(exp_dir * "/SWEpv.txt",SWEpv)
    # For functions that do not need to return a value (functions used only for some side effects), 
    # the Julia convention is to return the value nothing
    return nothing #SWEhat, GmeltHat, Ghat,  Phat, Ushat, G_pv, Gmelt_pv, U_pv, SWEpv #
    # return  # same as above (return keyword implicitly returns nothing, so it can be used alone)
end  # function
