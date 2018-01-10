
using GroupedErrors
using JLD,DataFrames,AxisArrays
using PlotUtils, RecipesBase, StatPlots
using LaTeXStrings, Measures
plotlyjs(gridcolor=:gray40,axiscolor=:gray50,
   textcolor=RGB(128/255,128/255,128/255),guidecolor=RGB(128/25,128/255,128/255),
    guidefontfamily="DejaVu Sans",guidefontsize=10,guidefontcolor=RGB(128/255,128/255,128/255) ,tickfont = ("DejaVu Sans",10))

try 
    mkdir("plots")
end

labbook_table = JLD.load("data/labbookTable.jld","df");

avg_data_dict = JLD.load("data/avgData.jld");

@load "data/statTables.jld"

cellPairs = sort(unique(labbook_table[isna.(labbook_table[:TAGS]),:cellToCell]))
genotypes = sort(unique(labbook_table[isna.(labbook_table[:TAGS]),:genotypeRegion]));

@recipe f{T}(aa::AxisArray{T,1}) =
(axes(aa,1)[:]-axes(aa,1)[findfirst(axes(aa,1)[:].>3.0)],aa)

@recipe f{T}(aa::AxisArray{T,3}) = 
(axes(aa,1)[:]-axes(aa,1)[findfirst(axes(aa,1)[:].>3.0)],
    reshape(aa,(size(aa,1)*size(aa,2),size(aa,3))))

@recipe f{T}(AA::Array{AxisArrays.AxisArray{T,1,Array{T,1},Tuple{AxisArrays.Axis{:time,StepRangeLen{Float64}}}},1}) =
          ([axes(aa,1)[:]-axes(aa,1)[findfirst(axes(aa,1)[:].>3.0)] for aa in AA],AA)

## Extract the traces and stats corresponding to a cell pair, number of pulses and range of runs (to select non drug runs)
function select_data(cellPair::Array{String},nPulses,data_dict,run_range,labbook,stats_per_run=stats_per_run)
    mini_labbook = labbook[findin(labbook[:cellToCell],cellPair),:]
    sort!(mini_labbook,cols=:keyEntry)
    pair_keys = mini_labbook[:keyEntry]
   
    mini_stats = stats_per_run[findin(stats_per_run[:experiment],pair_keys),:]
    if isa(run_range,Dict)
        mini_stats[:run_range] = [run_range[exp] for exp in mini_stats[:experiment]]
    else
        run_range =[run_range for i in 1:size(mini_stats,1)]
        mini_stats[:run_range] = run_range
    end

    mini_stats = 
    mini_stats[(in.(Array(mini_stats[:nPulses_median]),[nPulses])) .& [in(mini_stats[i,:runIdx],
                mini_stats[i,:run_range]) 
            for i in 1:size(mini_stats,1)],:]
  
    mini_labbook = mini_labbook[findin(pair_keys,mini_stats[:experiment]),:]
    sort!(mini_labbook,cols=(:keyEntry))
    genos=mini_labbook[:genotypeRegion]
    
    pair_keys = mini_stats[:experiment] 
    mini_data = [get(data_dict,pair_keys[i],0)[mini_stats[:runIdx][i]][1] for 
                                                         i in eachindex(pair_keys)];
    Dict("stats"=>mini_stats,"data"=>mini_data,"genotypes"=>genos)
end

function select_data(cellPair::String,nPulses,data_dict,run_range,labbook,stats_per_run=stats_per_run)
    select_data([cellPair],nPulses,data_dict,run_range,labbook,stats_per_run)
end

## Figure 2 of the paper
## Selected example pairs

pairs_for_figure2 = Dict([("vi","PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG1-8.s-EBt.b-DV_GA.b"),
                          ("v","PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl1.b-NO3PM.b"),
                          ("iv","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-FBl3.b-NO2D.b"),
                          ("iii","LAL.s-GAi.s-NO1i.b-to-EBIRP I-O-LAL.s"),
                          ("ii","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type1"),
                          ("i","PBG1-7.s-FBl2.s-LAL.b-cre.b-to-SMP.s-LAL.s-LAL.b.contra")
                         ])

function make_raw_plot(dataset;substract=false,scalebar=true,scaley=3,np=20,colorV=:experiment,traceW=3,kwargs...) ## substract is a keyword argument deciding if the
                                                          ## baseline fluorescence is to be substracted for display
    p = plot(textcolor=RGB(128/255,128/255,128/255),guidecolor=RGB(128/255,128/255,128/255);kwargs...)
    topvalue = 0
    colorsVIdx = unique(dataset["stats"][colorV])
    for i in 1:length(dataset["data"])
        d = copy(dataset["data"][i])
        if substract 
            d[:] = d .- dataset["stats"][:baseline_median][i]
        end
        plot!(p,d,color=find(colorsVIdx.==dataset["stats"][colorV][i])[1],legend=false,alpha=0.8,
            lw=traceW)
        topvalue = max(maximum(d),topvalue)
    end
    plot!(p,[0;0.033*np],[topvalue+1;topvalue+1],color=:gray80,fill=:gray80,alpha=0.4,fillrange=0,line=:path)
    if scalebar
        plot!(p,[6;8],[scaley-1;scaley-1],color=:gray50,line=:path,lw=3)
    end
    p
end

excPlot = make_raw_plot(select_data(pairs_for_figure2["i"],20,avg_data_dict,1:6,labbook_table),substract=true,
    xlims=(-2,10),framestyle=:axes,axis=:y)

inhPlot = make_raw_plot(select_data(pairs_for_figure2["ii"],20,avg_data_dict,1:6,labbook_table),scalebar=false,
    framestyle=:axes,axis=:y,xlims=(-2,10))

weakPlot = make_raw_plot(select_data(pairs_for_figure2["iii"],20,avg_data_dict,1:6,labbook_table),scaley=1.5,
    substract=true,scalebar=false,xlims=(-2,10),framestyle=:axes,axis=:y)

confPlot = make_raw_plot(select_data(pairs_for_figure2["iv"],20,avg_data_dict,1:6,labbook_table),
    substract=true,scalebar=false,framestyle=:axes,xlims=(-2,10),axis=:y)

reboundPlot = make_raw_plot(select_data(pairs_for_figure2["v"],20,avg_data_dict,1:6,labbook_table),substract=true,
   scalebar=false,framestyle=:axes,axis=:y,xlims=(-2,10))

nothingPlot = make_raw_plot(select_data(pairs_for_figure2["vi"],20,avg_data_dict,1:6,labbook_table),substract=true,
    scalebar=false,framestyle=:axes,axis=:y,xlims=(2,10))

direct_pair = "PBG1-7.s-FBl2.s-LAL.b-cre.b-to-PBG1-7.s-FBl2.s-LAL.b-cre.b"
direct = make_raw_plot(select_data(direct_pair,20,avg_data_dict,1:6,labbook_table),substract=true,
    framestyle=:axes,xlims=(-2,10))

figureResponseA = plot(excPlot,inhPlot,nothingPlot,weakPlot,confPlot,reboundPlot,layout = (2,3),
        title = ["i" "ii" "iii" "i" "ii" "iii"],titleloc=:left,top_margin=Measures.Length(:mm,30.0)
    ,size = (1000,600),bottom_margin=Measures.Length(:mm,15.0),link=:x)
## The top margin is to fit neuron schematics
savefig(figureResponseA,"plots/figureResponseAPlots.svg")

## Plotting the control experiments
null_pairs = convert(Array{String},unique(stats_per_run[stats_per_run[:expType].=="Non overlapping",:cellPair]))
figureResponseDi = make_raw_plot(select_data(null_pairs,20,avg_data_dict,1:6,labbook_table),substract=true,scalebar=true,framestyle=:axes,axis=:y,xlims=(-2,10),colorV=:cellPair,traceW=1,size=(400,350))
savefig(figureResponseDi,"plots/figureResponseDi.svg")

fullNamesDF = DataFrame(integNorm = "Normalized integral", 
                        integNormScaled = "Scaled normalized integral",
                        peakNorm = "Normalized peak",
                        peakFluo = "Peak fluorescence",
                        integ = "Integral",
                        decay_time = "Decay half time",
                        peakTime = "Rise time",
                        repeats_corr = "Repeat to repeat correlations",
                        between_runs_corr = "Correlations between experiments",
                        state_dependence_integral = "Integral to baseline correlation",
                        state_dependence = "Distance to baseline correlation",
                        #dose_slope_peak_norm= "Stimulation to normalized peak slope", 
                        #dose_slope_peak_normScaled= "Stimulation to normalized scaled peak slope", 
                        #dose_peak_norm = "Norm peak to stimulation correlation",  
                        distance = "Distance",
                        #baselineMAD = "Estimated baseline spread",
                        responding = "Responding");
 
makeStatHist = function(statDf,statN)
        stephist(statDf[statN],group=statDf[:expType],fill=true,alpha=0.4,normalize=:none,
            nbins=linspace(minimum(statDf[statN]),maximum(statDf[statN]),30),
            xlab=fullNamesDF[statN][1],
            color=["cornflowerblue" "coral" "green"])
end

makeStatHist! = function(statDf,pl,stat;kwargs...)
    stephist!(pl,statDf[stat],group=statDf[:expType],fill=true,alpha=0.4,normalize=:none,
            nbins=linspace(minimum(statDf[stat]),maximum(statDf[stat]),30),
            xlab="",
            color=["cornflowerblue" "coral" "green"];kwargs...)
end

##stats_per_pair_20 = stats_per_pair[stats_per_pair[:nPulses_median].==20,:]

statLayout = grid(2,2,widths=[0.9,0.1],heights=[0.1,0.9])
statLayout[1,2].attr[:blank]=true

statPlot = plot(layout=statLayout,size=(500,450),xlabel="",ylabel="",grid=true)

    makeStatHist!(stats_per_pair_20,statPlot,:repeats_corr,subplot=1,legend=:none,grid=false,link=:x,axis=:x,ticks=nothing)
makeStatHist!(stats_per_pair_20,statPlot,:integNormScaled,subplot=3,orientation=:h,legend=:none,grid=false,ylim=(-1.05,1.05),
    xlim = (0,18),axis=:y,ticks=nothing)

 scatter!(statPlot,stats_per_pair_20[:repeats_corr],stats_per_pair_20[:integNormScaled],group=stats_per_pair_20[:expType],
         msw=stats_per_pair_20[:signif5],msa=1,malpha=[0.8 0.3 0.3],msc=:gray30,
         ylab="Scaled normalized integral",xlab="Within-flies correlation",
      hover=stats_per_pair_20[:cellPair],subplot=2,msize=6,ylim=(-1.05,1.05),link=:x,color=["cornflowerblue" "coral" "green"])
    #scatter!([null_mean[2]],[null_mean[1]],markershape=:cross,subplot=2)
    

## There's still a hack in here, because orientation switching does weird things with plotlyjs (axis limits and 
## axis names don't match)

savefig(statPlot,"plots/figureResponsesDii.svg")
    
statHists = [makeStatHist(stats_per_pair_20,s) for s in names(fullNamesDF)[1:10]]

statsHistsGridBig = plot(statHists...,layout=(5,2),size=(800,1200),legend=:none)

savefig(statsHistsGridBig,"plots/statistics_histograms_SI.svg")

function remap(vals)
    remapped = zeros(Int, length(vals))
    for (i,x) in enumerate(uniqueTypesUsed)
        remapped[find(v -> v == x, vals)] = i
    end
    remapped
end

function getMat(I,J,column)
    mat = full(sparse(I, J, stats_per_pair_20[column],length(uniqueTypesUsed),length(uniqueTypesUsed)))
    mat[mat.==0.0] = NaN
    mat
end

Is = remap(stats_per_pair_20[:preNeuron])
Js = remap(stats_per_pair_20[:postNeuron])
statsMatrices = Dict(
"distanceN" => getMat(Is,Js,:distanceNorm)
)

matGuesses = transpose(full(sparse(Is, Js, stats_per_pair_20[:expType].=="Overlapping",
        length(uniqueTypesUsed),length(uniqueTypesUsed))))

function makeMatrixPlot(statUsed)
    mat = statsMatrices[statUsed]
    matVals = filter(x->!isnan(x),mat)
    grad = ColorGradient(ColorGradient(:bluesreds).colors,[0,-minimum(matVals)/(maximum(matVals)-minimum(matVals)),1.0])

    themat = heatmap(uniqueTypesUsed, uniqueTypesUsed, mat,aspect_ratio=1,
         color=grad,xrotation=90,size=(1500,1500),yguide="Presynaptic candidate",xguide="Postsynaptic candidate",
         xticks = (0.5:1:length(uniqueTypesUsed)-0.5,uniqueTypesUsed),
         yticks = (0.5:1:length(uniqueTypesUsed)-0.5,uniqueTypesUsed),
         )
    plot!(themat,identity,line=(:grey80,:dot),lab="",legend=false)
    scatter!(themat,collect(findn(matGuesses))-0.5...,m=(:black),msw=0,malpha=0.5,lab="")
    themat
end

matDistance = makeMatrixPlot("distanceN")

savefig(matDistance,"plots/matDistance.svg")

    linesToType = readtable("LinesAndTypes.csv")

shortPre = [linesToType[linesToType[:Type_Description].==n,:New_Type_Name][1] for n in stats_per_run[:preNeuron]]
shortPost = [linesToType[linesToType[:Type_Description].==n,:New_Type_Name][1] for n in stats_per_run[:postNeuron]]
stats_per_run[:shortPairName] = shortPre .* " to " .* shortPost
    
    
baselineDists = @df stats_per_run[stats_per_run[:preDrug],:] boxplot(:shortPairName,:baseline_median,
    size=(1500,400),ylims=(0,10),group=:expType,ylabel="Single run baseline",
    color=["cornflowerblue" "coral" "green"],whisker_width=0.5,xticks=length(:preDrug),
    linecolor=:gray50,markersize=2,alpha=0.8,malpha=0.6,xrotation=45,xtickfont = ("DejaVu Sans",6))

baselineDistsSummary = @df stats_per_run[stats_per_run[:preDrug],:] violin(:expType,:baseline_median,
    size=(600,400),ylabel="Single run baseline",ylims=(0,10),legend=:none,color="gray80")

stateDependenceSummary = @df stats_per_pair_20 boxplot(:globalSignif,
    :state_dependence_integral,
    size=(600,400),
    ylabel="Baseline to integral correlation",
    #group= :expType,
    #label=["Non overlapping" "Overlapping" "Self stimulation"],
    whisker_width=0.5,
    bar_width=0.7,
    #xrotation=20,
    fillcolor=:gray80,
    linecolor=:gray50,
    xticks=([-1,0,1],["Inhibition","Non significant","Excitation"]),
    xlim=(-1.6,1.5),
    legend=:none)

colorsF = Dict("i"=>"coral","ii"=>"coral","iii"=>"coral","iv"=>"green")
ylimsF = Dict("i"=>(0,5),"ii"=>(-1,0.2),"iii"=>(0,0.2),"iv"=>(0,6))
baselineCorr = [stats_per_pair_20[(stats_per_pair_20[:cellPair].==x) .& (
                           stats_per_pair_20[:nPulses_median].==20),:state_dependence] 
           for x in [[pairs_for_figure2[i] for i in ["i","ii","iii"]]...,direct_pair]]
figure2B = [@df select_data(pairs_for_figure2[x],20,avg_data_dict,1:6,labbook_table,
            stats_per_run
            )["stats"] scatter(:baseline_median,:integral_to_peak_median ,ylab="",
             mcolor = colorsF[x],
            xlab="",xlims=(0,),ylims=ylimsF[x]) for x in ["i","ii","iii"]]

push!(figure2B,
      @df select_data(direct_pair,20,avg_data_dict,1:6,
            labbook_table,stats_per_run)["stats"] scatter(:baseline_median,
        :integral_to_peak_median ,xlab="",ylab="",mcolor = colorsF["iv"]
        ,xlims=(0,1.7)))

plot(figure2B...,layout=(1,4),size=(1500,400),legend=false,msw=0,
    ylab=["Integral to peak" "" "" ""],
    xlab="Baseline value",
    title=["i" "ii" "iii" "iv"],titleloc=:left)

l=@layout [a
           b c
           d e f g]
baselineSIFig = plot(baselineDists,baselineDistsSummary,stateDependenceSummary,figure2B...,layout=l,size=(1500,1300),
    margin=Measures.Length(:mm,10.0),title=["A" "B" "C" "Di" "ii" "iii" "iv"],titleloc=:left,legend=[:right :none :right :none :none :none])

savefig(baselineSIFig,"plots/baselineSIFig.svg")

doseRespPlot = @> stats_per_pair begin 
    @splitby (_.signif20,_.cellPair)
    @x _.nPulses_median
    @y _.integNormScaled
    @set_attr :color _[1] == 1 ? :red : _[1]==0 ? :gray70 : :blue
    @plot plot(legend=:none,linealpha=0.3,size=(600,800),ylabel="Scaled normalized response integral",xlabel="Number of stimulation pulses",xticks=[1;5;10;20;30],xlims=(0,31))
end

@> stats_per_pair begin 
    @splitby _.signif20
    @across _.nPulses_median
    @x _.nPulses_median
    @y _.integNormScaled (mean,sem)
    @plot plot!(doseRespPlot,lw=4,color=[:blue :gray70 :red])
end

savefig(doseRespPlot,"plots/doseRespSI.svg")

#Blink.AtomShell.install()

@load "data/drugTables.jld"
interpData = JLD.load("data/interpolatedData.jld")

drugStats[:globalSignif]=stats_per_pair_20[:globalSignif][[findin(stats_per_pair_20[:cellPair],[s])[1] for s in drugStats[:cellPair]]];

function makePairDrugPlots(df,cp,findFunc=findmin)
    preData = []
    drugData = []
    postData = []
    pairPlot = plot(layout=(1,2),size=(500,500),legend=:none,
                    right_margin=[5mm 5mm 5mm 30mm 5mm 5mm 5mm 5mm])
    for expe in unique(df[df[:cellPair].==cp,:experiment])
        subDf = df[df[:experiment].==expe,:]
      
        runPre = subDf[subDf[:timeToDrug].<=2,:runIdx]
        runDrug = subDf[(10.<subDf[:timeToDrug].<15),:runIdx]
        runWash = subDf[(size(subDf,1)-1):size(subDf,1),:runIdx]
        # Within fly means
        push!(preData,squeeze(mean(cat(2,[interpData[expe][i][1] for i in runPre]...),2),2))
        push!(drugData,squeeze(mean(cat(2,[interpData[expe][i][1] for i in runDrug]...),2),2))
        push!(postData,squeeze(mean(cat(2,[interpData[expe][i][1] for i in runWash]...),2),2))
    end
    np = 30
    # Between flies means
    preDataM = squeeze(mean(cat(2,preData...),2),2)
    drugDataM = squeeze(mean(cat(2,drugData...),2),2)
    postDataM = squeeze(mean(cat(2,postData...),2),2)

    preDataS = squeeze(std(cat(2,preData...),2)./length(preData),2)
    drugDataS = squeeze(std(cat(2,drugData...),2)./length(preData),2)
    postDataS = squeeze(std(cat(2,postData...),2)./length(preData),2)

    plot!(preDataM,ribbon=preDataS,subplot=2,color=:cornflowerblue,xlabel="Time (s)")
    plot!(drugDataM,ribbon=drugDataS,subplot=2,color=:coral)
    plot!(postDataM,ribbon=postDataS,subplot=2,color=:gray70)
    
    topvalue = maximum([maximum(preDataM.+preDataS),maximum(postDataM.+postDataS),maximum(drugDataM.+drugDataS)])+1
    plot!([0;0.033*np],[topvalue;topvalue],color=:gray80,fill=:gray80,alpha=0.4,
        fillrange=0,line=:path,subplot=2)
    
    @> df[df[:cellPair].==cp,:] begin
    @splitby _.experiment
    @x _.timeToDrug
    @y _.integral_to_peak_scaled #_.distanceNorm 
    @set_attr :title cp 
    @plot scatter!(legend=:none,
        size=(600,800),ylabel="Response integral",xlabel="Time to drug",subplot=1,msw=0)
end
    pairPlot
end

#mecaEffectPlots = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[:cellPair]))]

mecaISP = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[mecadf[:preNeuron].=="PBG2-9.b-IB.s.SPS.s" ,:cellPair]))]

mecaColu = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[(mecadf[:preNeuron].=="PBG2-9.s-EBt.b-NO1.b.Type1").| (mecadf[:preNeuron].=="PBG1-7.s-FBl2.s-LAL.b-cre.b") .| (mecadf[:preNeuron].=="PBG1-8.b-EBw.s-DV_GA.b") .| (mecadf[:preNeuron].=="PBG1-8.s-EBt.b-DV_GA.b"),:cellPair]))]

mecaOthers = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[(mecadf[:preNeuron].=="EB.w-AMP.d-D_GAsurround").| (mecadf[:preNeuron].=="PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b") ,:cellPair]))]

mecaPlotsColu = plot(mecaColu...,layout=(5,2),size=(1000,1200),margin=3mm,
      bottom_margin=15mm,titlefont=font("DejaVu Sans",10),
    titlecolor=RGB(128/255,128/255,128/255),legend=:none)

savefig(mecaPlotsColu,"plots/mecaColuSI.svg")

mecaPlotsISP = plot(mecaISP...,plot(),layout=(2,2),size=(1000,480),margin=3mm,
      bottom_margin=15mm,titlefont=font("DejaVu Sans",10),
    titlecolor=RGB(128/255,128/255,128/255),legend=:none)

savefig(mecaPlotsColu,"plots/mecaISPSI.svg")

mecaPlotsOther = plot(mecaOthers...,layout=(1,2),size=(1000,240),margin=3mm,
      bottom_margin=15mm,titlefont=font("DejaVu Sans",10),
    titlecolor=RGB(128/255,128/255,128/255),legend=:none)

savefig(mecaPlotsOther,"plots/mecaOtherSI.svg")

#picroEffectPlots = [makePairDrugPlots(picrodf,cp,findmax) for cp in sort(unique(picrodf[:cellPair]))]

picroInhib = [makePairDrugPlots(picrodf,cp,findmax) for cp in ["EBIRP I-O-LAL.s-to-PBG1-8.b-EBw.s-DV_GA.b","EBORP O-I-GA-Bulb-to-PBG1-8.b-EBw.s-DV_GA.b","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type1","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type2","PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl3.b-NO2V.b"]]

picroControl = [makePairDrugPlots(picrodf,cp,findmax) for cp in ["PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-EBt.b-NO1.b.Type1","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-FBl1.b-NO3PM.b","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.b-IB.s.SPS.s","PBG2-9.b-IB.s.SPS.s-to-PBG2-9.s-EBt.b-NO1.b.Type1","PBG2-9.b-IB.s.SPS.s-to-PBG1-8.b-EBw.s-DV_GA.b"    ,"PBG2-9.b-IB.s.SPS.s-to-PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b" ,"PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-EBt.b-NO1.b.Type2","EB.w-AMP.d-D_GAsurround-to-PBG1-8.b-EBw.s-DV_GA.b"]]

picroPlotsInhib = plot(picroInhib...,plot([]),layout=(3,2),size=(1000,720),margin=3mm,
      bottom_margin=15mm,titlefont=font("DejaVu Sans",10),titlecolor=RGB(128/255,128/255,128/255),legend=:none)

savefig(picroPlotsInhib,"plots/picroInhibSI.svg")

picroPlotsControl = plot(picroControl...,layout=(4,2),size=(1000,960),margin=3mm,
      bottom_margin=15mm,titlefont=font("DejaVu Sans",10),titlecolor=RGB(128/255,128/255,128/255),legend=:none)

savefig(picroPlotsControl,"plots/picroControlSI.svg")

inhibPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl3.b-NO2V.b"
excitPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-EBt.b-NO1.b.Type1"
mixedPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG1-8.b-EBw.s-DV_GA.b"

deltaInhib = make_raw_plot(select_data(inhibPair,20,avg_data_dict,1:6,labbook_table),scalebar=false,substract=false,
    substract=false,xlims=(-2,10),xlabel = "Time (seconds)", ylabel = "Fluorescence", title = "A",titleloc=:left)

deltaExcit = make_raw_plot(select_data(excitPair,20,avg_data_dict,1:6,labbook_table),
    scalebar=false,substract=false,xlims=(-2,10),xlabel = "Time (seconds)", 
    ylabel = "Fluorescence", title = "B",titleloc=:left)

deltaMixed = [make_raw_plot(select_data(mixedPair,p,avg_data_dict,1:6,labbook_table),
        scalebar=false,substract=false,np=p,xlims=(-2,10)) for p in [5,10,20,30]]

deltaMixedP = plot(deltaMixed...,layout=(1,4),size=(1500,500),title=["C i" "ii" "iii" "iv"],titleloc=:left,
    ylabel=["Fluorescence" "" "" ""],xlabel=["" "" "" "Time (seconds)"],margin=Measures.Length(:mm,10.0))

deltaL = @layout [g h
                  z]
deltaFig = plot(deltaInhib,deltaExcit,deltaMixedP,layout=deltaL,size=(1500,1000),
    bottom_margin=Measures.Length(:mm,20.0))

savefig(deltaFig,"delta7SI.svg")
