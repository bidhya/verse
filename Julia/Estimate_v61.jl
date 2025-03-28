using Random
using JuMP
using Ipopt
using DelimitedFiles
Random.seed!(1234)  # seed for reproducibility

function blender(i, j, SWEprior, Pprior, Gprior, SCFinst, AirT, logDir, exp_dir, twindow)
  """
  inputs
  ==============
  i,j are pixel locations
  SWEprior, Pprior, Gprior are prior estimates of SWE [m], Precipitation [m] and surface heat flux [W/m2]
  SCFinst is measured snow covered fraction
  AirT is air temperature [K]
  logdir and exp_dir are directories where log and data are written
    
  variable sizes
  =============== 
  SWE and SCF are storage terms, so will be length nt. 
  P, Melt, and G are flux terms, so will be length nt-1
  for input, all variables are length nt, so just use Pprior[1:nt-1]. Note, Gprior currently unused
  for output in section 4, will output fluxes as nt, with the last element set to 0.    

  """
    # 0 handle variable sizes
    nt = length(SCFinst)
    Pprior = Pprior[1:nt-1]
    
    # 1 Smooth SCF observations    
    # twindow = 5
    SCFobs = smoothdata(SCFinst,twindow,nt,"mean");    
    # twindow = 60
    SCF_smooth_season = smoothdata(SCFinst,60,nt,"mean");

    # 2 Define hyperparameters
    tmelt,tmelt_smooth,SWEmax,SWEmin_global,Meltmax,σP,σSWE,k,Melt0,L=define_hyperparameters(SCF_smooth_season,nt,Pprior,SWEprior,AirT, SCFobs)

    # 3 Solve
    m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>5000))
    # set_silent(m)
    # define variables and bounds
    @variable(m, SWEmin_global <= SWE[i=1:nt] <= SWEmax[i], start=SWEprior[i]);
    @variable(m, Precip[i=1:nt-1]>=0. ,start=Pprior[i]);
    @variable(m, 0. <= Melt[i=1:nt-1]<= Meltmax);
    @variable(m, Mcost[i=1:nt-1] >=0);
    # define constraints
    for i in 1:nt-1
      @NLconstraint(m,Mcost[i]==L/(1+exp( -k*(Melt[i]-Melt0)   ) ))
    end
    for i in 1:nt-1
      @constraint(m,SWE[i+1]==SWE[i]+Precip[i]-Melt[i])
    end
    # define objective function
    @objective(m,Min,sum((Precip-Pprior).^2 ./σP.^2) + sum((SWE-SWEprior ).^2 ./ σSWE.^2) + sum(Mcost.^2));
    log_file =  "$logDir/Pix_$(i)_$(j)_$(twindow).txt"  # original
    # solve
    redirect_stdio(stdout=log_file, stderr=log_file) do
      optimize!(m)
    end
    # print(termination_status(m))
    
    # 4 extract
    NODATAvalue = -9999
    SWEhat=JuMP.value.(SWE);
    Phat=zeros(nt,1)
    Phat[1:nt-1] = JuMP.value.(Precip);           
    Melt_hat = zeros(nt,1)
    Melt_hat[1:nt-1] = JuMP.value.(Melt);
    
    Δt = 86400; #s/day
    ρw = 1000; #density of water
    Lf = 0.334E6; #Latent heat of fusion J/kg    
    GmeltHat = Melt_hat/Δt*Lf*ρw
    
    Ghat = ones(nt,1)*NODATAvalue;    
    Ushat = ones(nt,1)*NODATAvalue;
    G_pv = ones(nt,1)*NODATAvalue;
    U_pv = ones(nt,1)*NODATAvalue;
    Gmelt_pv = ones(nt,1)*NODATAvalue;
    SWEpv = ones(nt,1)*NODATAvalue;
    
    # 5 output
    out_vars = hcat(SWEhat, GmeltHat, Ghat, Phat, Ushat, G_pv, Gmelt_pv, U_pv, SWEpv, SCFobs)
    writedlm("$(exp_dir)/Pix_$(i)_$(j)_$(twindow).txt", out_vars)
    
    # 6 clean up
    m = nothing
    GC.gc()
    return nothing    
end

function smoothdata(SCF_inst,twindow,nt,smoothfunc)
    SCF_smooth=zeros(nt,1);
    for i=1:nt
        istart = trunc(Int,i-round(twindow/2))
        iend = trunc(Int,i+round(twindow/2))        
        # if i < twindow || i > nt-twindow
        if istart < 1 || iend > nt
            SCF_smooth[i]=0
        else            
            if smoothfunc == "mean"
                SCF_smooth[i] = mean(SCF_inst[istart:iend])
            elseif smoothfunc=="median"
                SCF_smooth[i] = median(SCF_inst[istart:iend])
            end
        end
    end
    
    return SCF_smooth
end

function define_uncertainty(Pprior,SWEprior,AirT,SCFobs,nt,tmelt_smooth)
    # convert air temperature K-> C
    AirT=AirT.-273.15
    
    # 2.2.1 Precipitation Uncertainty
    RelPUnc = 0.3; #[-] this applies to cumulative precipitation
    # Uncertainty for accumulation . precip is size nt-1
    σP = zeros(nt-1,1)
    σPmin = 0.001
    Pthresh = 0.001
    for i=1:nt-1
      if Pprior[i]<Pthresh
        σP[i]=σPmin
      else
        σP[i]=Pprior[i]*RelPUnc
      end
    end    
    # adjust uncertainty to apply to the number of snow days
    nsnowday = 0
    Tprecip_thresh=1.5
    for i=1:nt-1
        if Pprior[i]>Pthresh && AirT[i] < Tprecip_thresh
            nsnowday += 1
        end
    end
    if nsnowday > 0
        σP = σP*sqrt(nsnowday);    
    end
    
    # 2.2.2 SWE Uncertainty
    fSWE = 0.4
    σSWE = SWEprior*fSWE;
    σSWEmin = 0.01
    σSWEmax = 10
    for i=1:nt
        # if SWEprior[i]>0 && SCFobs[i]==0
        if SCFobs[i]==0 || tmelt_smooth[i]>.1    
            σSWE[i] = σSWEmax
        end
        if σSWE[i] < σSWEmin
            σSWE[i] = σSWEmin
        end
    end   
    
    return σP, σSWE
end

function define_hyperparameters(SCF_smooth_season,nt,Pprior,SWEprior,AirT, SCFobs)
    """
    list of hyperparameters
    =======================
    twindow for smoothing SCF for snow on/off constraint - defined in main
    twindow for smoothing SCF for identifying melt times - defined in main
    ΔSCFthresh
    twindow for smoothing melt times
    SWEmax_global
    SWEmin_global
    Meltmax    
    k,Melt0,L
    σPmin - defined in define_uncertainty
    Pthresh -  defined in define_uncertainty
    Tprecip_thresh -  defined in define_uncertainty
    fSWE -  defined in define_uncertainty
    σSWEmin -  defined in define_uncertainty
    σSWEmax -  defined in define_uncertainty        
    
    """
    
    # 2.1 Define prior estimates
    # 2.1.1 Define times when snow is melting 
    tmelt = zeros(nt,1)
    ΔSCFthresh = -0.01
    for i=2:nt
        if SCF_smooth_season[i]-SCF_smooth_season[i-1] < ΔSCFthresh
            tmelt[i] = 1
        end
    end    
    # twindow = 30
    tmelt_smooth = smoothdata(tmelt,30,nt,"mean");
    
    # 2.2 Extreme / limit values
    # 2.2.1 SWE
    SWEmax_global = 5;
    SWEmin_global = 1.0e-6; #1/1000 mm
    # define SWEmax as a function of time and of SCF
    #    set SWEmax to 0 if SCF is low
    SWEmax = zeros(nt,1)
    for i = 1:nt
      if SCFobs[i] == 0
        SWEmax[i] = SWEmin_global
      else
        SWEmax[i] = SWEmax_global
      end
    end    
    #2.2.2 Melt
    Meltmax = 0.075;
    
    # 2.3 Define uncertainty
    σP,σSWE = define_uncertainty(Pprior,SWEprior,AirT,SCFobs,nt,tmelt_smooth)
    # 2.4 Melt cost function parameters
    k = 500
    Melt0 = 0.05
    L = 1
    
    return tmelt,tmelt_smooth,SWEmax,SWEmin_global,Meltmax,σP,σSWE,k,Melt0,L
end