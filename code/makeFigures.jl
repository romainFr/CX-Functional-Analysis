
using GroupedErrors
using JLD2,DataFrames,AxisArrays,CSV,FileIO
using PlotUtils, RecipesBase, StatPlots, StatsBase
using LaTeXStrings, Measures
theme(:wong,gridcolor=:gray40,axiscolor=:gray50,
   textcolor=RGB(128/255,128/255,128/255),guidecolor=RGB(128/25,128/255,128/255),
    guidefontfamily="DejaVu Sans",guidefontsize=10,guidefontcolor=RGB(128/255,128/255,128/255) ,tickfont = ("DejaVu Sans",10),legendfont=("DejaVu Sans",10,RGB(128/255,128/255,128/255)),foreground_color=RGB(128/255,128/255,128/255),background_color=RGBA(1.0,1.0,1.0,0.0))
plotlyjs()

include("code/functions/plot-utilities.jl")

mkpath("plots")

@load "data/labbookTable.jld2" labbook

avg_data_dict = load("data/avgData.jld2")

@load "data/statTables.jld2"

@load "data/drugTables.jld2"
interpData = load("data/interpolatedData.jld2")

cellPairs = sort(unique(labbook[ismissing.(labbook[:TAGS]),:cellToCell]))
genotypes = sort(unique(labbook[ismissing.(labbook[:TAGS]),:genotypeRegion]));

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
        title = ["A i" "ii" "iii" "B i" "ii" "iii"],titleloc=:left,top_margin=30mm
    ,size = (800,500),bottom_margin=7mm,link=:x)
## The top margin is to fit neuron schematics
savefig(figureResponseA,"plots/figureResponseAPlots.svg")
PlotlyJS.savefig(figureResponseA.o,"plots/figureResponseAPlots.html",js=:remote)

## Plotting the control experiments
null_pairs = convert(Array{String},unique(stats_per_run[stats_per_run[:expType].=="Non overlapping",:cellPair]))

figureResponseDLayout = @layout [i grid(2,2,widths=[0.9,0.1],heights=[0.1,0.9])]
figureResponseDLayout[1,2][1,2].attr[:blank]=true
figureResponseD = plot(layout=figureResponseDLayout,size=(800,400),xlabel="",ylabel="",legend=false,titlepos=:left)

make_raw_plot!(figureResponseD,select_data(null_pairs,20,interpData,1:6,labbook),substract=true,scalebar=true,framestyle=:axes,axis=:y,colorV=:cellPair,traceW=1,right_margin=10mm,top_margin=10mm,label="",subplot=1,title="D i")
makeStatHist!(stats_per_pair_20,figureResponseD,:between_runs_corr,subplot=2,grid=false,link=:x,axis=:x,ticks=nothing,top_margin=15mm,label="",title="ii")

makeStatHist!(stats_per_pair_20,figureResponseD,:integNormScaled,subplot=4,orientation=:h,grid=false,axis=:y,ticks=nothing,ylim=(-1.05,1.05),xlim = (0,18),legend=(0.67,1.02))
 
@df stats_per_pair_20 scatter!(figureResponseD,:between_runs_corr,:integNormScaled,
                               group=(:signif1,:expType),
                               color=[1 2 1 2 3],
                               msw=[0 0 1 1 1],msa=1,
                               malpha=0.4,msc=:gray30,
                               ylab="Scaled normalized integral",xlab="Between-flies correlation",
                               hover=:cellPair,subplot=3,msize=6,ylim=(-1.05,1.05),link=:x
                               ,label=["" "" "" "" "Significant response"])

## There's still a hack in here, because orientation switching does weird things with plotlyjs (axis limits and 
## axis names don't match)

savefig(figureResponseD,"plots/figureResponsesD.svg")
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
 

statHists = [makeStatHist(stats_per_pair_20,names(fullNamesDF)[1]),[makeStatHist(stats_per_pair_20,s,label="") for s in names(fullNamesDF)[2:10]]...]

statsHistsGridBig = plot(statHists...,layout=(5,2),size=(700,700),legend=(0.8,0.97),margin=5mm)

savefig(statsHistsGridBig,"plots/statistics_histograms_SI.svg")
savefig(statsHistsGridBig,"plots/statistics_histograms_SI.pdf")
PlotlyJS.savefig(statsHistsGridBig.o,"plots/statistics_histograms_SI.html",js=:remote)


Is = remap(stats_per_pair_20[:preNeuron])
Js = remap(stats_per_pair_20[:postNeuron])
statsMatrices = Dict(
"distanceN" => getMat(Is,Js,:distanceNorm)
)

matGuesses = transpose(full(sparse(Is, Js, stats_per_pair_20[:expType].=="Overlapping",
        length(uniqueTypesUsed),length(uniqueTypesUsed))))

matDistance = 
makeMatrixPlot("distanceN",size=(1000,1000),title="B",title_location=:left,top_margin=5mm,tickfontsize=9,legend=(0.1,0.2))

savefig(matDistance,"plots/matDistance.svg")
PlotlyJS.savefig(matDistance.o,"plots/matDistance.html",js=:remote)

## Baseline effects
baselineDists = @df stats_per_run[stats_per_run[:preDrug],:] boxplot(:cellPair,:baseline_median,
    size=(800,300),ylims=(0,10),group=:expType,ylabel="Single run baseline",
                                                                     whisker_width=0.5,xticks=[],
                                                                     linecolor=:gray50,markersize=2,alpha=0.8,malpha=0.6,
                                                                     hover=:cellPair,legend=(0.35,0.92))

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

colorsF = Dict("i"=>2,"ii"=>2,"iii"=>2,"iv"=>3)
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

savefig(baselineSIFig,"plots/baselineSIFig.svg")
savefig(baselineSIFig,"plots/baselineSIFig.pdf")
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

savefig(doseRespPlot,"plots/doseRespSI.svg")
savefig(doseRespPlot,"plots/doseRespSI.pdf")
PlotlyJS.savefig(doseRespPlot.o,"plots/doseRespSI.html",js=:remote)

#Blink.AtomShell.install()



#drugStats[:globalSignif]=stats_per_pair_20[:globalSignif][[findin(stats_per_pair_20[:cellPair],[s])[1] for s in drugStats[:cellPair]]];

#mecaEffectPlots = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[:cellPair]))]

mecaISP = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[mecadf[:preNeuron].=="PBG2-9.b-IB.s.SPS.s" ,:cellPair]))]

mecaColu = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[(mecadf[:preNeuron].=="PBG2-9.s-EBt.b-NO1.b.Type1").| (mecadf[:preNeuron].=="PBG1-7.s-FBl2.s-LAL.b-cre.b") .| (mecadf[:preNeuron].=="PBG1-8.b-EBw.s-DV_GA.b") .| (mecadf[:preNeuron].=="PBG1-8.s-EBt.b-DV_GA.b"),:cellPair]))]

mecaOthers = [makePairDrugPlots(mecadf,cp) for cp in sort(unique(mecadf[(mecadf[:preNeuron].=="EB.w-AMP.d-D_GAsurround").| (mecadf[:preNeuron].=="PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b") ,:cellPair]))]

mecaOthers[1].subplots[2].series_list[1][:label] = "Control"
mecaOthers[1].subplots[2].series_list[2][:label] = "Mecamylamine"
mecaOthers[1].subplots[2].series_list[3][:label] = "Wash"

mecaLayout = @layout [a{0.01h}
                      grid(5,2){0.635h}
                      b{0.01h}
                      grid(2,2){0.265h}
                      c{.01h}
                      grid(1,2){0.07h}]

mecaPlots = plot(plot(title="A",title_location=:left,framestyle=:none),
                 mecaColu...,
                 plot(title="B",title_location=:left,framestyle=:none),
                 mecaISP...,plot(),
                 plot(title="C",title_location=:left,framestyle=:none),
                 mecaOthers...,
                 layout=mecaLayout,size=(800,1500),legend=(0.75,0.2))

savefig(mecaPlots,"plots/mecaPlots.svg")
PlotlyJS.savefig(mecaPlots.o,"plots/mecaPlots.html",js=:remote)


picroInhib = [makePairDrugPlots(picrodf,cp) for cp in ["EBIRP I-O-LAL.s-to-PBG1-8.b-EBw.s-DV_GA.b","EBORP O-I-GA-Bulb-to-PBG1-8.b-EBw.s-DV_GA.b","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type1","LAL.s-GAi.s-NO1i.b-to-PBG2-9.s-EBt.b-NO1.b.Type2","PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl3.b-NO2V.b"]]

picroControl = [makePairDrugPlots(picrodf,cp) for cp in ["PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-EBt.b-NO1.b.Type1","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.s-FBl1.b-NO3PM.b","PBG1-8.b-EBw.s-DV_GA.b-to-PBG2-9.b-IB.s.SPS.s","PBG2-9.b-IB.s.SPS.s-to-PBG2-9.s-EBt.b-NO1.b.Type1","PBG2-9.b-IB.s.SPS.s-to-PBG1-8.b-EBw.s-DV_GA.b"    ,"PBG2-9.b-IB.s.SPS.s-to-PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b" ,"PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-EBt.b-NO1.b.Type2","EB.w-AMP.d-D_GAsurround-to-PBG1-8.b-EBw.s-DV_GA.b"]]

picroInhib[1].subplots[2].series_list[1][:label] = "Control"
picroInhib[1].subplots[2].series_list[2][:label] = "Picrotoxin"
picroInhib[1].subplots[2].series_list[3][:label] = "Wash"

picroLayout = @layout [a{0.01h}
                       grid(3,2){0.42h}
                       b{0.01h}
                       grid(4,2){0.56h}]

picroPlots = plot(plot(title="A",title_location=:left,framestyle=:none),
                  picroInhib...,plot(),
                  plot(title="B",title_location=:left,framestyle=:none),
                  picroControl...,
                  layout=picroLayout,size=(1000,2100),legend=(0.75,0.6)
                  )

savefig(picroPlots,"plots/picroPlots.svg")
PlotlyJS.savefig(picroPlots.o,"plots/picroPlots.html",js=:remote)

inhibPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-FBl3.b-NO2V.b"
excitPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG2-9.s-EBt.b-NO1.b.Type1"
mixedPair = "PB18.s-GxΔ7Gy.b-PB18.s-9i1i8c.b-to-PBG1-8.b-EBw.s-DV_GA.b"

deltaInhib = make_raw_plot(select_data(inhibPair,20,
                                       #avg_data_dict
                                       interpData,1:6,labbook),scalebar=false,substract=false,
    substract=false,xlabel = "Time (seconds)", ylabel = "Fluorescence", title = "A",titleloc=:left,right_margin=3mm)
 
deltaExcit = make_raw_plot(select_data(excitPair,20,interpData,1:6,labbook),
    scalebar=false,substract=false,xlabel = "Time (seconds)", 
    ylabel = "Fluorescence", title = "B",titleloc=:left,right_margin=3mm)

deltaMixed = [make_raw_plot(select_data(mixedPair,p,interpData,1:6,labbook),
        scalebar=false,substract=false,np=p) for p in [5,10,20,30]]

deltaMixedP = plot(deltaMixed...,layout=(1,4),size=(1500,500),title=["C i" "ii" "iii" "iv"],titleloc=:left,
    ylabel=["Fluorescence" "" "" ""],xlabel=["" "" "" "Time (seconds)"],right_margin=3mm,left_margin=3mm)

deltaL = @layout [g h
                  z]
deltaFig = plot(deltaInhib,deltaExcit,deltaMixedP,layout=deltaL,size=(800,800),
    bottom_margin=7mm,top_margin=5mm)

savefig(deltaFig,"plots/delta7SI.svg")
PlotlyJS.savefig(deltaFig.o,"plots/delta7SI.html",js=:remote)
