
using GroupedErrors
using JLD2,DataFrames,AxisArrays,CSV,FileIO
using PlotUtils, RecipesBase, StatPlots, StatsBase
using LaTeXStrings, Measures, CSV, Rsvg

myPalette = [colorant"#94a6fd",
             colorant"#841ea4",
             colorant"#69b7c5",
             colorant"#085782",
             colorant"#da73f8",
             colorant"#21a645",
             colorant"#285d28",
             colorant"#9db46a"]
            
#theme()
plotlyjs(palette=PlotThemes.expand_palette(colorant"white", myPalette; lchoices=linspace(57,57,1), cchoices=linspace(100,100,1)),gridcolor=:gray40,axiscolor=:gray50,
   textcolor=RGB(128/255,128/255,128/255),guidecolor=RGB(128/25,128/255,128/255),
    guidefontfamily="DejaVu Sans",guidefontsize=10,guidefontcolor=RGB(128/255,128/255,128/255) ,tickfont = ("DejaVu Sans",10),legendfont=("DejaVu Sans",10,RGB(128/255,128/255,128/255)),foreground_color=RGB(128/255,128/255,128/255),background_color=RGBA(1.0,1.0,1.0,0.0))

include("code/functions/plot-utilities.jl")

mkpath("plots")

@load "data/labbookTable.jld2" labbook

avg_data_dict = load("data/avgData.jld2")

@load "data/statTables.jld2"

@load "data/drugTables.jld2"
interpData = load("data/interpolatedData.jld2")

linesToType = CSV.read("LinesAndTypes.csv",strings=:raw)

cellPairs = sort(unique(labbook[ismissing.(labbook[:TAGS]),:cellToCell]))
genotypes = sort(unique(labbook[ismissing.(labbook[:TAGS]),:genotypeRegion]));

mecadf[:shortPair] = 
linesToType[[findfirst(linesToType[Symbol("Type Description")],pn) for pn in mecadf[:preNeuron]],Symbol("New Type Name")].*" to ".*linesToType[[findfirst(linesToType[Symbol("Type Description")],pn) for pn in mecadf[:postNeuron]],Symbol("New Type Name")]

picrodf[:shortPair] = linesToType[[findfirst(linesToType[Symbol("Type Description")],pn) for pn in picrodf[:preNeuron]],Symbol("New Type Name")].*" to ".*linesToType[[findfirst(linesToType[Symbol("Type Description")],pn) for pn in picrodf[:postNeuron]],Symbol("New Type Name")]

stats_per_run[:shortPair] = linesToType[[findfirst(linesToType[Symbol("Type Description")],pn) for pn in stats_per_run[:preNeuron]],Symbol("New Type Name")].*" to ".*linesToType[[findfirst(linesToType[Symbol("Type Description")],pn) for pn in stats_per_run[:postNeuron]],Symbol("New Type Name")]

uniqueTypesShort = linesToType[[findfirst(linesToType[Symbol("Type Description")],pn) for pn in uniqueTypesUsed],Symbol("New Type Name")]
uniqueTypesShort = convert(Array{String,1},uniqueTypesShort)
## Figure "Responses" of the paper
## Selected example pairs

pairs_for_figureResp = Dict([("vi","PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG1-8.s-EBt.b-DV_GA.b"),
                          ("v","PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl1.b-NO3PM.b"),
                          ("iv","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-FBl3.b-NO2D.b"),
                          ("iii","LAL.s-GAi.s-NO1i.b-to-EBIRP I-O-LAL.s"),
                          ("ii","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type1"),
                          ("i","PBG1-7.s-FBl2.s-LAL.b-cre.b-to-SMP.s-LAL.s-LAL.b.contra")
                         ])


excPlot = make_raw_plot(select_data(pairs_for_figureResp["i"],20,interpData,1:6,labbook),substract=true,framestyle=:axes,axis=:y,link=:x)

inhPlot = make_raw_plot(select_data(pairs_for_figureResp["ii"],20,interpData,1:6,labbook),scalebar=false,
    framestyle=:axes,axis=:y,link=:x)

weakPlot = make_raw_plot(select_data(pairs_for_figureResp["iii"],20,interpData,1:6,labbook),scaley=1.5,
    substract=true,scalebar=false,framestyle=:axes,axis=:y,link=:x)

confPlot = make_raw_plot(select_data(pairs_for_figureResp["iv"],20,interpData,1:6,labbook),
    substract=true,scalebar=false,framestyle=:axes,axis=:y,link=:x)
 
reboundPlot = make_raw_plot(select_data(pairs_for_figureResp["v"],20,interpData,1:6,labbook),substract=true,
   scalebar=false,framestyle=:axes,axis=:y,link=:x,yticks=[0;1.5;3])

nothingPlot = make_raw_plot(select_data(pairs_for_figureResp["vi"],20,interpData,1:6,labbook),substract=true,
    scalebar=false,framestyle=:axes,axis=:y,link=:x)
   
direct_pair = "PBG1-7.s-FBl2.s-LAL.b-cre.b-to-PBG1-7.s-FBl2.s-LAL.b-cre.b"
direct = make_raw_plot(select_data(direct_pair,20,interpData,1:6,labbook),substract=true,
    framestyle=:axes)

figureResponseA = plot(excPlot,inhPlot,nothingPlot,weakPlot,confPlot,reboundPlot,layout = (2,3),
        title = ["A i" "ii" "iii" "B i" "ii" "iii"],ylabel=["Fluorescence (ΔF/F₀)" "" "" "" "" "" ""],titleloc=:left,top_margin=30mm
    ,size = (800,500),bottom_margin=7mm,link=:x)
## The top margin is to fit neuron schematics
PlotlyJS.savefig(figureResponseA.o,"plots/figureResponseAPlots.svg")
PlotlyJS.savefig(figureResponseA.o,"plots/figureResponseAPlots.html",js=:remote)

## Plotting the control experiments
null_pairs = convert(Array{String},unique(stats_per_run[stats_per_run[:expType].=="Non overlapping",:cellPair]))
null_pairs_short = convert(Array{String},unique(stats_per_run[stats_per_run[:expType].=="Non overlapping",:shortPair]))
nullPlots = [make_raw_plot(select_data(np,20,interpData,1:6,labbook),substract=true,scalebar=false,framestyle=:axes,colorV=:cellPair,right_margin=5mm,top_margin=5mm,label="",lw=3) for np in null_pairs]

nullPlots[10].subplots[1].attr[:yaxis][:guide] = "Fluorescence (ΔF/F₀)"

nullPlot = plot(nullPlots...,layout=(7,3),size=(750,1400),link=:x,title=reshape([null_pairs_short...,"","",""],1,length(null_pairs_short)+3),xlab="Time (s)")

PlotlyJS.savefig(nullPlot.o,"plots/SIFigureNonOverlapping.svg")
PlotlyJS.savefig(nullPlot.o,"plots/SIFigureNonOverlapping.html",js=:remote)

##
figureResponseDLayout = grid(2,2,widths=[0.9,0.1],heights=[0.1,0.9])
figureResponseDLayout[1,2].attr[:blank]=true
figureResponseD = plot(layout=figureResponseDLayout,size=(500,600),xlabel="",ylabel="",legend=false,titlepos=:left)

makeStatHist!(stats_per_pair_20,figureResponseD,:between_runs_corr,subplot=1,grid=false,link=:x,axis=:x,ticks=nothing,top_margin=15mm,label="",title="D",color=[1 6 9])

makeStatHist!(stats_per_pair_20,figureResponseD,:integNormScaled,subplot=3,orientation=:h,grid=false,axis=:y,ticks=nothing,ylim=(-1.05,1.05),xlim = (0,18),legend=(0.4,1.02),color=[1 6 9])
 
@df stats_per_pair_20 scatter!(figureResponseD,:between_runs_corr,:integNormScaled,
                               group=(:signif1,:expType),
                               color=[1 6 1 6 9],
                               msw=[0 0 1 1 1],msa=1,
                               malpha=0.4,msc=:gray30,
                               ylab="Scaled normalized integral",xlab="Between-flies correlation",
                               hover=:cellPair,subplot=2,msize=6,ylim=(-1.05,1.05),link=:x
                               ,label=["" "" "" "" "Significant response"])

## There's still a hack in here, because orientation switching does weird things with plotlyjs (axis limits and 
## axis names don't match)

PlotlyJS.savefig(figureResponseD.o,"plots/figureResponsesD.svg")
PlotlyJS.savefig(figureResponseD.o,"plots/figureResponsesD.html",js=:remote)



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
                        distance = "Distance",
                        responding = "Responding");
 

statHists = [makeStatHist(stats_per_pair_20,names(fullNamesDF)[1],color=[1 6 9]),[makeStatHist(stats_per_pair_20,s,label="",color=[1 6 9]) for s in names(fullNamesDF)[2:10]]...]

statsHistsGridBig = plot(statHists...,layout=(5,2),size=(700,700),legend=(0.8,0.97),margin=5mm)

PlotlyJS.savefig(statsHistsGridBig.o,"plots/statistics_histograms_SI.svg")
PlotlyJS.savefig(statsHistsGridBig.o,"plots/statistics_histograms_SI.html",js=:remote)


Is = remap(stats_per_pair_20[:preNeuron])
Js = remap(stats_per_pair_20[:postNeuron])
statsMatrices = Dict(
"distanceN" => getMat(Is,Js,:distanceNorm)
)

matGuesses = transpose(full(sparse(Is, Js, stats_per_pair_20[:expType].=="Overlapping",
        length(uniqueTypesUsed),length(uniqueTypesUsed))))

matDistance = 
       makeMatrixPlot("distanceN",uniqueTypesShort,size=(1000,1000),title_location=:left,top_margin=5mm,tickfontsize=9,legend=(0.02,-0.01))

PlotlyJS.savefig(matDistance.o,"plots/matDistance.svg")
PlotlyJS.savefig(matDistance.o,"plots/matDistance.html",js=:remote)

## Baseline effects
baselineDists = @df stats_per_run[stats_per_run[:preDrug],:] boxplot(:cellPair,:baseline_median,
    size=(800,300),ylims=(0,10),group=:expType,ylabel="Single run baseline",
                                                                     whisker_width=0.5,xticks=[],
                                                                     linecolor=:gray50,markersize=2,alpha=0.8,malpha=0.6,
                                                                     hover=:cellPair,legend=(0.35,0.92),color=[1 6 9])

baselineDistsSummary = @df stats_per_run[stats_per_run[:preDrug],:] violin(:expType,:baseline_median,
    size=(350,300),ylabel="Single run baseline",ylims=(0,10),legend=:none,color="gray80",label="")

stateDependenceSummary = @df stats_per_pair_20 boxplot(:globalSignif,
    :state_dependence_integral,
    size=(350,300),
    ylabel="Baseline to integral correlation",
    whisker_width=0.5,
    bar_width=0.7,
    fillcolor=:gray80,
    linecolor=:gray50,
    xticks=([-1,0,1],["Inhibition","Non significant","Excitation"]),
    xlim=(-1.6,1.5),
    legend=:none,label="")

colorsF = Dict("i"=>6,"ii"=>6,"iii"=>6,"iv"=>9)
ylimsF = Dict("i"=>(0,5),"ii"=>(-1,0.2),"iii"=>(0,0.2),"iv"=>(0,6))

figureExampleState = [@df select_data(pairs_for_figureResp[x],20,avg_data_dict,1:6,labbook,
            stats_per_run
            )["stats"] scatter(:baseline_median,:integral_to_peak_median ,ylab="",
             mcolor = colorsF[x],
            xlab="",xlims=(0,),ylims=ylimsF[x],label="") for x in ["i","ii","iii"]]

push!(figureExampleState,
      @df select_data(direct_pair,20,avg_data_dict,1:6,
            labbook,stats_per_run)["stats"] scatter(:baseline_median,
        :integral_to_peak_median ,xlab="",ylab="",mcolor = colorsF["iv"]
                                                    ,xlims=(0,1.7),label=""))

plot(figureExampleState...,layout=(1,4),legend=false,msw=0,
    ylab=["Integral to peak" "" "" ""],
    xlab="Baseline value",
    title=["i" "ii" "iii" "iv"],titleloc=:left,xrotation=45)

l=@layout [a
           b c
           d e f g]
baselineSIFig = plot(baselineDists,baselineDistsSummary,stateDependenceSummary,figureExampleState...,layout=l,size=(800,900),
    top_margin=5mm,bottom_margin=7mm,title=["A" "B" "C" "Di" "ii" "iii" "iv"],titleloc=:left,legend=(0.35,0.97))

PlotlyJS.savefig(baselineSIFig.o,"plots/baselineSIFig.svg")
PlotlyJS.savefig(baselineSIFig.o,"plots/baselineSIFig.html",js=:remote)

doseRespPlot = @> stats_per_pair begin 
    @splitby (_.signif20,_.cellPair)
    @x _.nPulses_median
    @y _.integNormScaled
    @set_attr :color _[1] == 1 ? :darkred : _[1]==0 ? :gray70 : :darkblue
    @set_attr :hover _[2]
    @plot plot(label="",linealpha=0.3,line=:line,size=(400,500),ylabel="Scaled normalized response integral",xlabel="Number of stimulation pulses",xticks=[1;5;10;20;30],xlims=(0,31))
end

@> stats_per_pair begin 
    @splitby _.signif20
    @across _.nPulses_median
    @x _.nPulses_median
    @y _.integNormScaled (mean,sem)
    @plot plot!(doseRespPlot,lw=4,color=[:darkblue :gray70 :darkred],legend=(0.15,0.97),label=["Excitatory responses" "No responses" "Inhibitory responses"])
end

PlotlyJS.savefig(doseRespPlot.o,"plots/doseRespSI.svg")
PlotlyJS.savefig(doseRespPlot.o,"plots/doseRespSI.html",js=:remote)

#Blink.AtomShell.install()



#drugStats[:globalSignif]=stats_per_pair_20[:globalSignif][[findin(stats_per_pair_20[:cellPair],[s])[1] for s in drugStats[:cellPair]]];


mecaISP = vcat([makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[mecadf[:preNeuron].=="PBG2-9.b-IB.s.SPS.s" ,:cellPair]))]...)

mecaColu = vcat([makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[(mecadf[:preNeuron].=="PBG2-9.s-EBt.b-NO1.b.Type1").| (mecadf[:preNeuron].=="PBG1-7.s-FBl2.s-LAL.b-cre.b") .| (mecadf[:preNeuron].=="PBG1-8.b-EBw.s-DV_GA.b") .| (mecadf[:preNeuron].=="PBG1-8.s-EBt.b-DV_GA.b"),:cellPair]))]...)

mecaOthers = vcat([makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[(mecadf[:preNeuron].=="EB.w-AMP.d-D_GAsurround").| (mecadf[:preNeuron].=="PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b") ,:cellPair]))]...)

mecaOthers[2].subplots[1].series_list[1][:label] = "Control"
mecaOthers[2].subplots[1].series_list[2][:label] = "Mecamylamine"
mecaOthers[2].subplots[1].series_list[3][:label] = "Wash"
mecaColu[17].subplots[1].attr[:yaxis][:guide] = "Scaled normalized integral"

mecaPlots = plot(
                 mecaColu...,
                 mecaISP...,plot(),plot(),
                 mecaOthers...,
                 layout=grid(8,4),size=(800,1700),legend=(0.75,0.15),link=:x,xlab=["Time to drug (min)" "Time (s)"])

PlotlyJS.savefig(mecaPlots.o,"plots/mecaPlots.svg")
PlotlyJS.savefig(mecaPlots.o,"plots/mecaPlots.html",js=:remote)


picroInhib = vcat([makePairDrugPlots(picrodf,cp) for cp in ["EBIRP I-O-LAL.s-to-PBG1-8.b-EBw.s-DV_GA.b","EBORP O-I-GA-Bulb-to-PBG1-8.b-EBw.s-DV_GA.b","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type1","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type2","PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl3.b-NO2V.b"]]...)

picroControl = vcat([makePairDrugPlots(picrodf,cp) for cp in ["PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-EBt.b-NO1.b.Type1","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-FBl1.b-NO3PM.b","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.b-IB.s.SPS.s","PBG2-9.b-IB.s.SPS.s-to-PBG2-9.s-EBt.b-NO1.b.Type1","PBG2-9.b-IB.s.SPS.s-to-PBG1-8.b-EBw.s-DV_GA.b"    ,"PBG2-9.b-IB.s.SPS.s-to-PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b" ,"PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-EBt.b-NO1.b.Type2","EB.w-AMP.d-D_GAsurround-to-PBG1-8.b-EBw.s-DV_GA.b"]]...)

picroInhib[2].subplots[1].series_list[1][:label] = "Control"
picroInhib[2].subplots[1].series_list[2][:label] = "Picrotoxin"
picroInhib[2].subplots[1].series_list[3][:label] = "Wash"
picroControl[1].subplots[1].attr[:yaxis][:guide] = "Scaled normalized integral"

picroPlots = plot(
                  picroInhib...,plot(),plot(),
                  picroControl...,
                  layout=grid(7,4),size=(800,1500),legend=(0.75,0.6),link=:x,xlab=["Time to drug (min)" "Time (s)"]
                  )

PlotlyJS.savefig(picroPlots.o,"plots/picroPlots.svg")
PlotlyJS.savefig(picroPlots.o,"plots/picroPlots.html",js=:remote)


## SI Figure : mixed effects of Delta7
inhibPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl3.b-NO2V.b"
excitPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-EBt.b-NO1.b.Type1"
mixedPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG1-8.b-EBw.s-DV_GA.b"

deltaInhib = make_raw_plot(select_data(inhibPair,20,
                                       #avg_data_dict
                                       interpData,1:6,labbook),scalebar=false,substract=false,
    substract=false,xlabel = "Time (s)", ylabel = "Fluorescence (ΔF/F₀)", title = "A",titleloc=:left,right_margin=3mm)
 
deltaExcit = make_raw_plot(select_data(excitPair,20,interpData,1:6,labbook),
    scalebar=false,substract=false,xlabel = "Time (s)", 
    ylabel = "Fluorescence (ΔF/F₀)", title = "B",titleloc=:left,right_margin=3mm)

deltaMixed = [make_raw_plot(select_data(mixedPair,p,interpData,1:6,labbook),
        scalebar=false,substract=false,np=p) for p in [5,10,20,30]]

deltaMixedP = plot(deltaMixed...,layout=(1,4),size=(1500,500),title=["C i" "ii" "iii" "iv"],titleloc=:left,
    ylabel=["Fluorescence (ΔF/F₀)" "" "" ""],xlabel=["" "" "" "Time (s)"],right_margin=3mm,left_margin=3mm)

deltaL = @layout [g h
                  z]
deltaFig = plot(deltaInhib,deltaExcit,deltaMixedP,layout=deltaL,size=(800,800),
    bottom_margin=7mm,top_margin=30mm)

PlotlyJS.savefig(deltaFig.o,"plots/delta7SI.svg")
PlotlyJS.savefig(deltaFig.o,"plots/delta7SI.html",js=:remote)

## SI figure : effect of genotypes
## Extracting the pairs where this happens
nGenoDf = by(labbook,:cellToCell,df -> length(unique(df[:LexA])))
genoEffectPairs = nGenoDf[nGenoDf[:x1].>1,:cellToCell]

genoEffectPlots = [make_raw_plot(select_data(cp,20,interpData,1:6,labbook),scalebar=false,colorV=:genotype,legendV=true,substract=false,pairTitle=true) for cp in genoEffectPairs]
genoEffectPlots[2].subplots[1].attr[:yaxis][:guide] = "Fluorescence (ΔF/F₀)"

genoEffectPlot = plot(genoEffectPlots[1],plot(),genoEffectPlots[2:end]...,layout=(3,2),size=(800,1000),link=:x,xlab="Time (s)")
PlotlyJS.savefig(genoEffectPlot.o,"plots/genotypeEffectSI.svg")


## SI figure : effect of region imaged
nRegionDf = by(labbook,:genotype,df -> DataFrame(cellToCell =df[:cellToCell][1],nReg = length(unique(df[:Region]))*(size(df,1)>3))) 
regionEffectPairs = nRegionDf[nRegionDf[:nReg].>1,:cellToCell]

regionEffectPlots = [make_raw_plot(select_data(cp,20,interpData,1:6,labbook),scalebar=false,colorV=:region,legendV=true,substract=false,pairTitle=true) for cp in regionEffectPairs]
regionEffectPlots[10].subplots[1].attr[:yaxis][:guide] = "Fluorescence (ΔF/F₀)"

regionEffectPlot = plot(regionEffectPlots...,layout=(7,3),size=(1600,3500),link=:x,xlab="Time (s)")
PlotlyJS.savefig(regionEffectPlot.o,"plots/regionEffectSI.svg")
