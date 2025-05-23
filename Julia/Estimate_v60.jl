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
# v55 Fixing error due to missing days of data due to Polar nights
# v58 New updates by jack (Nov 2023)
# Jan 29, 2024: Due to error on obj function σWRFG becoming zero, fixed the minimum σWRFG to 25 
# Jun 21, 2024: No change for 1km run here. MODIS with UINT8 did not work likely because of missing due to no-data background and polar nights.
# Nov 5, 2024: simpler heat flux parameterization

using Random
using JuMP
using Ipopt
using DelimitedFiles
Random.seed!(1234)  # seed for reproducibility

function fix_modis(SCF)
    # Define length of array (365)
    nt = length(SCF)

    # %% 1.1 Calculate deltaSCF
    # We calculate the value ΔSCF which is defined as
    # Σ[ abs[SCF(i)-SCF(i+1)] + abs[SCF(i+1)-SCF(i+2)] ]
    # Depending on value of ΔSCF we do one of three things
    
    for i in 150:nt-2  # Python uses 0-based indexing, so 150 in MATLAB is 149 in Python
        deltaSCF = abs(SCF[i] - SCF[i+1]) + abs(SCF[i+1] - SCF[i+2])
        # print(deltaSCF)
        if deltaSCF == 2  # If ΔSCF == 2, that means the SCF went from 1→0→1 which is bad
            tmp = SCF[i:i+2]  #.copy()
            tmp[tmp .== 0] .= 0.667  # In this case we find the 0 SCF day (which should be in the middle) and then we set that days SCF = 0.667
            SCF[i:i+2] = tmp
            # println(2)
        elseif deltaSCF > 1.5  # Elif ΔSCF > 1.5, we set any 0 SCF day to 1/3 of value Σ[abs(ΔSCF(i:1+2))]
            tmp = SCF[i:i+2]  #.copy()
            tmp[tmp .== 0] .= deltaSCF/3
            SCF[i:i+2] = tmp
            println(1.5)
        elseif deltaSCF > 0.9  # Elif ΔSCF > 0.9, we set any 0 SCF day to 1/2 of value Σ[abs(ΔSCF(i:1+2))]
            tmp = SCF[i:i+2] #.copy()
            tmp[tmp .== 0] .= deltaSCF/2
            SCF[i:i+2] = tmp
            # print(0.9)
        end
    end
    # %% 1.2 Find final day of snow off
    # After the loop above finishes, find the final day of snow off

    # extra to check for 
    stopIdx = 0
    for i in 150:nt
        if SCF[i] == 0  # Find a day with zero SCF
            global stopIdx = i  # Set counter value to whatever idx i is
            numSCF = sum(SCF[i:nt])  # Sum all remaining SCF values
            if numSCF == 0  # If there is no more SCF for the rest of the year, break
                break
            end
        end
    end
    # Add a single step down day to the SCF timeseries
    # This ensures if the last day of SCF is above 50% snow cover
    # We add a single extra day to the timeseries where we cut that value down
    # by 50% to add an easier downramp for the melt timeseries
    # use try/catch to avoid error if stopIdx was not assigned above
    if stopIdx > 1 && SCF[stopIdx-1] > 0.49  # error if stopIdx was not assigned above; hence, wrapping in try block
        SCF[stopIdx] = 0.5 * SCF[stopIdx-1]
    end
    return SCF
end

function sigmaG(opt, WRFSWE, MSCF, Gmelt_pv, nt, fG=0.3, σWRFGmin=1)
    """
    if we pass WRFSWE and MSCF then opt is not really required.

    This function returns the sigma_G and Gmelt_prior
        fg: fraction of G
    Inputs:
    =======
        WRFSWE: LIS data
        MSCF: MODIS data
    """
    # fG = 0.3  # fraction of G
    Gmelt_prior=zeros(nt,1)
    σWRFG=zeros(nt,1) #.+ 1 # default of 1 everywhere
    if opt == 1 # default
        for i=1:nt
            if i>240 && i<280
                Gmelt_prior[i]=50
                σWRFG[i]=50
            else
                Gmelt_prior[i]=0
                σWRFG[i]=5
            end
        end
    elseif opt == 2
        # TODO: What will be threshold for no snow cover? == 0. or some small value?
        for i=1:nt
            if WRFSWE[i]>0. && MSCF[i] >0. # both LIS and MODIS says snowy
                Gmelt_prior[i] = Gmelt_pv[i]
                σWRFG[i] = fG * abs(Gmelt_pv[i])
            elseif WRFSWE[i]<=0. && MSCF[i]<=0. # both LIS and MODIS says not snowy
                Gmelt_prior[i] = 0
                σWRFG[i] = 1
            elseif WRFSWE[i] <= 0. && MSCF[i] <= 0. # LIS not snowy but MODIS snowy
                if i > 120  # Late season; approximated here around start of February
                    Gmelt_prior[i] = 50
                    σWRFG[i] = fG * abs(Gmelt_pv[i])
                else
                    Gmelt_prior[i] = 5
                    σWRFG[i] = fG * abs(Gmelt_pv[i])
                end
            elseif WRFSWE[i] > 0. && MSCF[i] <= 0. # LIS snowy but MODIS says not snow
                Gmelt_prior[i] = 0
                σWRFG[i] = 1
            end
            if σWRFG[i] == 0  # to guard against zero σWRFG error when dividing by it in objective function
                # Due to division by σWRFG when can become zero; propogated from Gmelt_pv being 0.0
                σWRFG[i] = 1
            end
        end
    end
    # set minimum limit on σWRFG
    for i=1:nt
	if σWRFG[i] < σWRFGmin
	    σWRFG[i]=σWRFGmin
        end
    end

	
    return Gmelt_prior, σWRFG
end

# using Dates, CSV, DataFrames
# function blender(DataDir, exp_dir, WRFSWE, WRFP, WRFG, MSCF, AirT)
# function blender(out_folder, i, j, WRFSWE, WRFP, WRFG, MSCF, AirT)
function blender(i, j, WRFSWE, WRFP, WRFG, MSCF, AirT, logDir, exp_dir, opt, σWRFGmin, fix_modis_flag=0)
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
    if fix_modis_flag==1
        MSCF = fix_modis(MSCF)  # apply MODIS fix developed by Jack (Feb 04, 2025)
    end

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
    Pmax=0.2 # 170 (rounded to 200) mm in one day of SWE corresponds to Sierra Nevada snowfall record
    ρnew=100 #density of new snow
    z0=0.01 # roughness length used in SDC
    mf=1.0 # melt factor used in SDC
    ρWRF=350 #note: previous versions used time-varying snow density but that is not currently exported

    # 1.1.1. read in σWRFG based on SCF

    #currently unused until further testing


    # 1.2 physical parameters
    ρw=1000 #density of water
    Lf=0.334E6 #Latent heat of fusion J/kg
    cs_ice=2102 #specific heat capacity for ice J/kg/K

    # 1.3 Fill in missing SCF [Feb 23, 2023]
    MissingSCFData=zeros(nt,1)
    for i=1:nt
        if ismissing(MSCF[i])
            MissingSCFData[i]=1
            MSCF[i]=1  # June 20, 2024: cannot convert a value to missing for assignment. When using uint8 based data.
        end
    end

    #= new heat flux uncertainty parameterization
    # New Updates fro Jack (Nov, 2023)
    σWRFG = zeros(nt,1) .+ 15  # maybe 15 is not being used, but we still need to initialze the array
    σWRFG_rel = 0.5
    for i=1:nt
        if MissingSCFData[i]==1 # missing data check
	        σWRFG[i] = 1e9
        elseif MSCF[i]>=0.1 && WRFSWE[i]>=0.1  # if both are snow covered
            σWRFG[i] = abs(WRFG[i]) * σWRFG_rel  # when WRFG == zero, σWRFG also becomes zero which creates problem of objective function
            # Jan 29, 2024 To solve this problem in obj function: Expression contains an invalid NaN constant. This could be produced by `Inf - Inf`.
            if σWRFG[i] < 25  # or == 0 because this is what was actually causing problem.  
                σWRFG[i] = 25  # σWRFG_rel  #abs(WRFG[i] + eps(Float32)) * σWRFG_rel
            end
        elseif MSCF[i]<0.1 && WRFSWE[i]<0.1 # if both not snowy
            σWRFG[i] = 25
        else
            σWRFG[i] = 500 #if they disagree, then don’t use prior in cost function
        end
    end
    =#

    # 1.4 Match up SWE and MSCF
    for i=1:nt
        if MSCF[i]==0  # Nov 04, 2022: "ERROR: LoadError: TypeError: non-boolean (Missing) used in boolean context" because some days there was no MODIS data
            WRFSWE[i] = 0
        end
    end
    
    # 2. compute useful variables and set up  arrays
    # 2.1 define SWEmax: upper limit of SWE, based on observed SCF
    SWEmax=zeros(nt,1)
    for i=1:nt
        if MSCF[i]==0  # TypeError: non-boolean (Missing) used in boolean context
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
    # Jan 30, 2023: Modify σWRFP based on number snowy days  
    nsnowday = 0  # number of snow days
    for i=1:nt
        if WRFP[i]>0.001 && AirT[i] < 1.5
            # there is precipitation and temperature is likely to convert it to snow.
            nsnowday += 1
        end
    end
    if nsnowday > 0
        # to prevent σWRFP from becoming zero if nsnowday=0.
        σWRFP = σWRFP .* sqrt(nsnowday * 0.5)
    end

    # change precipitation to snowfall, by setting precipitation to zero if air temp is too high
    for i=1:nt 
	if WRFP[i]>0 && AirT[i] > 1.5
	    WRFP[i]=0
	end
    end
    
    # # New-BNY For Arctic night [Feb 14, 2023] <-- Superceded by updates from Jack in Nov 2023.
    # # Set the vector to 15 everythere unless SCF is undefined for Arctic Nights, then set to large number
    # σWRFG=zeros(nt,1) .+ 15  # default of 15 everywhere
    # for i=1:nt
    #     # if MSCF[i] == 0.
    #     # if ismissing(MSCF[i])
    #     if MissingSCFData[i]==1
    #         # TODO: Check whether these values are NaN, missing, 0 or sth else
    #         σWRFG[i] = 1e9
    #     # else
    #     #     σWRFG[i] = 15
    #     end
    # end

    
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
        if Gmelt_pv[i]>0. && MSCF[i] >0.  # TypeError: non-boolean (Missing) used in boolean context
            #G_pv[i]=Gmelt_pv[i]/MSCF[i]
	    G_pv[i]=Gmelt_pv[i]
            U_pv[i]=0.
        else
            G_pv[i]=WRFG[i]
            #U_pv[i]=U_pv[i-1]+WRFP[i-1]*AirT[i-1]*ρw*cs_ice + G_pv[i-1]*MSCF[i-1]*Δt  # MethodError: Cannot `convert` an object of type Missing to an object of type Float64
	    U_pv[i]=U_pv[i-1]+WRFP[i-1]*AirT[i-1]*ρw*cs_ice + G_pv[i-1]*Δt  # 
            if U_pv[i]>0. || SWEpv[i]==0.
                U_pv[i]=0.
            end
        end
    end

    #=
    for i=2:nt
        if Gmelt_pv[i] > Gmax
            σWRFG[i]=1.0e9
        end
    end
    =#

    # this is a much simpler parameterization of prior heat flux and its uncertainty (Nov5, 2024)
    # TODO: Make this a function
    # TODO call the function here. we will have Gmelt_prior, σWRFG
    # opt = 1 #1 # 2 
    Gmelt_prior, σWRFG = sigmaG(opt, WRFSWE, MSCF, Gmelt_pv, nt, 0.5, σWRFGmin)  # 100 added σWRFGmin=1; run with σWRFGmin set to 10 and 20 and 100 …. keep fG=0.5 for all runs
    # println("Min Max: Gmelt_pv $(minimum(Gmelt_pv))  $(maximum(Gmelt_pv)) and σWRFG: $(minimum(σWRFG))   $(maximum(σWRFG))")

    # Gmelt_prior=zeros(nt,1)
    # σWRFG=zeros(nt,1)
    # for i=1:nt
    #     if i>240 && i<280
    #         Gmelt_prior[i]=50
    #         σWRFG[i]=50
    #     else
    #         Gmelt_prior[i]=0
    #         σWRFG[i]=5
    #     end
    # end
			

    # 3.2 Solve for the posterior using prior valid
    m = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>5000))
    #@variable(m, SWEmin <= SWE[i=1:nt] <= SWEmax[i],start=SWEpv[i] )
    @variable(m, SWEmin <= SWE[i=1:nt] <= SWEmax[i],start=WRFSWE[i] )
    @variable(m, 0 <= Precip[i=1:nt]<= Pmax,start=WRFP[i])
    #@variable(m,  G[i=1:nt] , start=G_pv[i])
    @variable(m,  G[i=1:nt] , start=WRFG[i])
    #@variable(m, 0 <= Gmelt[i=1:nt] <= Gmax, start=Gmelt_pv[i])
    @variable(m, 0 <= Gmelt[i=1:nt] <= Gmax, start=WRFG[i])	
    #@variable(m, Us[i=1:nt] <=0, start=U_pv[i])
    @variable(m, Us[i=1:nt] <=0, start=0.)
    #@objective(m,Min,sum( (SWE-SWEpv).^2 ./σWSWE.^2)+sum((Precip-WRFP).^2 ./σWRFP.^2)+
    #                   sum((G-G_pv).^2 ./σWRFG.^2) )
    #@objective(m,Min,sum((Precip-WRFP).^2 ./σWRFP.^2)+ sum((Gmelt-Gmelt_pv).^2 ./σWRFG.^2) )
    @objective(m,Min,sum((Precip-WRFP).^2 ./σWRFP.^2)+ sum((Gmelt-Gmelt_prior).^2 ./σWRFG.^2) )  #  Expression contains an invalid NaN constant. This could be produced by `Inf - Inf`.
    for i in 1:nt-1
        @constraint(m,SWE[i+1]==SWE[i]+Precip[i]-Gmelt[i]*Δt/Lf/ρw)
        #@NLconstraint(m,Us[i+1]==Us[i]+(1-(tanh(Us[i]/10000)+1))*G[i]*MSCF[i]*Δt+
        @NLconstraint(m,Us[i+1]==Us[i]+(1-(tanh(Us[i]/10000)+1))*G[i]*Δt+			
                                Precip[i]*AirT[i]*ρw*cs_ice) #m x K x kg/m3 x J/kg/K
    end
    @constraint(m,Us[1]==0) 
    @constraint(m,Us[nt]==0)
    for i in 1:nt
        #@NLconstraint(m,Gmelt[i]==G[i]*MSCF[i]*(tanh(Us[i]/10000)+1))
	@NLconstraint(m,Gmelt[i]==G[i]*(tanh(Us[i]/10000)+1))
    end
    
    # optimize!(m)
    # BY Create output folder; no error if the folder already exist
    # exp_dir = "$out_folder/outputs_txt"    # To save text outputs for each pixel
    # exp_dir = "$(tmp_txtDir)/Pix_$(i)_$(j)"
    # logDir = "logs"  #"$out_folder/logs"    # To save text outputs for each pixel
    # log_file =  "$logDir/Pix_$(i)_$(j).txt"  # original
    log_file =  "$logDir/Pix_$(i)_$(j)_$(σWRFGmin)_$(fix_modis_flag).txt"
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

    # out_vars = hcat(SWEhat, GmeltHat, Ghat, Phat, Ushat, G_pv, Gmelt_pv, U_pv, SWEpv)  # 9 columns original
    # writedlm("$(exp_dir)/Pix_$(i)_$(j).txt", out_vars)
    out_vars = hcat(SWEhat, GmeltHat, Ghat, Phat, Ushat, G_pv, Gmelt_pv, U_pv, SWEpv, Gmelt_prior, σWRFG)    # only for testing
    out_vars = hcat(out_vars, WRFSWE, WRFP, WRFG, MSCF, AirT)  # append inputs as well so we have one file for plotting and analysis
    writedlm("$(exp_dir)/Pix_$(i)_$(j)_$(σWRFGmin)_$(fix_modis_flag).txt", out_vars)  # Feb 05, 2025: _$(σWRFGmin)_$(fix_modis_flag) -> new suffix added to output file name for prototyping 

    # writedlm(exp_dir * "/SWE.txt",SWEhat)
    # writedlm(exp_dir * "/Gmelt.txt",GmeltHat)
    # writedlm(exp_dir * "/G.txt",Ghat)
    # writedlm(exp_dir * "/Precip.txt",Phat)
    # writedlm(exp_dir * "/Us.txt",Ushat)
    # writedlm(exp_dir * "/Gpv.txt",G_pv)
    # writedlm(exp_dir * "/Gmeltpv.txt",Gmelt_pv)
    # writedlm(exp_dir * "/Upv.txt",U_pv)
    # writedlm(exp_dir * "/SWEpv.txt",SWEpv)
    # To address growing memory issue but also causing: signal (11.1): Segmentation fault 
    m = nothing
    GC.gc()
    # For functions that do not need to return a value (functions used only for some side effects), 
    # the Julia convention is to return the value nothing
    return nothing
    # return SWEhat, GmeltHat, Ghat,  Phat, Ushat, G_pv, Gmelt_pv, U_pv, SWEpv #    
    # return  # same as above (return keyword implicitly returns nothing, so it can be used alone)
end  # function
