using JLD2,DataFrames,AxisArrays,FileIO,CSV
using Interpolations
using StatsBase,Distributions
using JSON
using DataStructures
using Distances,Bootstrap,StatsBase

include("functions/fluoRunUtilities.jl")
## Load the labbook, the lines description table and the fluorescence data
@load "data/labbookTable.jld2" labbook
linesToType = CSV.read("LinesAndTypes.csv")
full_data_dict = load("data/rawData.jld2");

## An per run average version of the fluorescence data (to make figures)
avg_data_dict = map(full_data_dict) do full  
    (k,full) = full
    meanSignals = [(AxisArray(vec(mean(f[1].data,3)),axes(f[1],Axis{:time})),f[2]) for 
        f in full]
    Pair(k,meanSignals)    
end;
save("data/avgData.jld2",avg_data_dict)

## Interpolating those averaged runs (useful to compute correlations)
interpolated_data_dict = map(avg_data_dict) do dct
        (k,val) = dct
        newAvg = [interpolate_run(run,3:0.1:5) for run in val]     ## The two seconds following the stimulation
        Pair(k,newAvg)
end;

large_interpolated_data_dict = map(avg_data_dict) do dct
        (k,val) = dct
        newAvg = [interpolate_run(run,1:0.1:12) for run in val]     ## The two seconds following the stimulation
        Pair(k,newAvg)
end;

save("data/interpolatedData.jld2",large_interpolated_data_dict)

keyEntries = labbook[:keyEntry];
fullStats = compileStats(full_data_dict,keyEntries)

## Calculate the medians of the statistics
stats_per_run = aggregate(fullStats,[:experiment,:runIdx],median)
stats_per_run[:nPulses_median] = convert(Array{Int},stats_per_run[:nPulses_median]);

## Add some metadata columns 
stats_per_run[:cellPair]=""
stats_per_run[:genotype]=""
stats_per_run[:preNeuron]=""
stats_per_run[:postNeuron]=""
stats_per_run[:preDrug]=true
stats_per_run[:Drug]=""
stats_per_run[:timeToDrug]= -Inf
for i in 1:size(stats_per_run,1)
    stats_per_run[i,:cellPair] = labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),:cellToCell]
    stats_per_run[i,:genotype] = labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),:genotypeRegion]
    stats_per_run[i,:preNeuron] = labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),:cellPre]
    stats_per_run[i,:postNeuron] = labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),:cellPost]
    if (!ismissing(labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),:Drug]))
        stats_per_run[i,:preDrug] = (labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),
            :timesToDrug][stats_per_run[i,:runIdx]])>Dates.Second(0)
        stats_per_run[i,:timeToDrug] = -Int64(Dates.value(labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),
            :timesToDrug][stats_per_run[i,:runIdx]]))/60000
        stats_per_run[i,:Drug]=labbook[findfirst(keyEntries.==stats_per_run[i,:experiment]),:Drug]
    end
end

## Calculate the number of runs per pair and select pairs with at least 3 runs
ns = by(stats_per_run,:cellPair,df -> DataFrame(n = length(unique(df[:experiment]))))
fullExps = ns[ns[:n].>2,:cellPair];

stats_per_run = stats_per_run[in.(Array(stats_per_run[:cellPair]),[fullExps]),:]

## List the cell types tested
uniqueTypesUsed = unique(vcat(stats_per_run[:preNeuron],stats_per_run[:postNeuron]))
sort!(uniqueTypesUsed,by= v->linesToType[findfirst(linesToType[Symbol("Type Description")].==v),:Supertype] )

## Parse the neuron types to define their pre/post regions
possibleNeuropiles = ["PB","FB","EB","NO","GA","LAL","rub","BU"];
neuronTypes = OrderedDict(uniqueTypesUsed[i] =>
                          Dict("innervates" => filter(x -> contains(uniqueTypesUsed[i],x),
                                                                possibleNeuropiles), 
                               "pre"=> split(linesToType[findfirst(linesToType[:,Symbol("Type Description")].==uniqueTypesUsed[i]),
                Symbol("Pre regions")],","),
                               "post" => split(linesToType[findfirst(linesToType[:,Symbol("Type Description")].==uniqueTypesUsed[i]),
                Symbol("Post regions")],","),
                               "pre_fine"=> split(linesToType[findfirst(linesToType[:,Symbol("Type Description")].==uniqueTypesUsed[i]),
                Symbol("Pre regions fine")],","),
                                "post_fine" => split(linesToType[findfirst(linesToType[:,Symbol("Type Description")].==uniqueTypesUsed[i]),
                Symbol("Post regions fine")],","),
                                "short_name" => linesToType[findfirst(linesToType[:,Symbol("Type Description")].==uniqueTypesUsed[i]),
            Symbol("New Type Name")]) for i in 1:length(uniqueTypesUsed))

## Using the annotation, establish if there's a potential overlap for every pair 
stats_per_run[:overlapping]=
[(length(intersect(neuronTypes[ty]["pre_fine"],
        neuronTypes[tyPost]["post_fine"]))>0) for (ty,tyPost) in zip(stats_per_run[:preNeuron],
                                                                     stats_per_run[:postNeuron])]

## Is the same neuron used for recording and stimulation ?
stats_per_run[:self]= (stats_per_run[:preNeuron].==stats_per_run[:postNeuron]);

stats_per_run[:expType] = "Non overlapping"
stats_per_run[stats_per_run[:overlapping],:expType] = "Overlapping"
stats_per_run[stats_per_run[:self],:expType] = "Self stimulation"
stats_per_run[:overlapping] =  stats_per_run[:expType].=="Overlapping"; ## Excluding self activation

stats_per_run[:integNorm_scaled] =  scaleResponse(stats_per_run[:integNorm_median],trimming=true)
stats_per_run[:integral_to_peak_scaled] =  scaleResponse(stats_per_run[:integral_to_peak_median],trimming=true)

## Drug experiments
# Selecting mecamylamine runs    
mecadf = stats_per_run[(stats_per_run[:Drug].=="Mecamylamine") .& (stats_per_run[:timeToDrug].>-5),:]

# Same thing for picrotoxin
picrodf = stats_per_run[(stats_per_run[:Drug].=="Picrotoxin") .& (stats_per_run[:timeToDrug].>-5),:]

## Calculate the stats for each pair
stats_per_pair = by(stats_per_run[stats_per_run[:preDrug],:],[:cellPair,:nPulses_median],category_stats)

## Scaled responses
stats_per_pair[:integNormScaled] = scaleResponse(stats_per_pair[:integNorm])
stats_to_use = [:integNormScaled,
                 :between_runs_corr
                 #:repeats_corr
               ];

addDistances!(stats_per_pair,stats_to_use,:integNormScaled)
## Add a signed significance to be used in summary diagrams/matrices
stats_per_pair[:globalSignif] =  sign.(stats_per_pair[:distance]).*stats_per_pair[:signif1]
stats_per_pair[stats_per_pair[:globalSignif].==-0.0,:globalSignif]=0
stats_per_pair[:distanceNorm] = stats_per_pair[:distance]./maximum(trim(abs.(stats_per_pair[stats_per_pair[:nPulses_median].<=20,:distance]),
                                                                        prop=0.01));

## Using the 20 pulses significance results for everything here
stats_per_pair_20 = stats_per_pair[stats_per_pair[:nPulses_median].==20,:]
stats_per_pair[:signif20] = stats_per_pair_20[indexin(stats_per_pair[:cellPair],
                                                      stats_per_pair_20[:cellPair]),:globalSignif]

#stats_per_pair_20 = join(stats_per_pair_20,doseRespDF,on=:cellPair,kind=:left)

## Export the stats (to be used by the plotting scripts)
@save "data/statTables.jld2" stats_per_run stats_per_pair uniqueTypesUsed stats_per_pair_20

@save "data/drugTables.jld2" mecadf picrodf

## Export all the tables to csv
CSV.write("data/stats_per_run.csv",stats_per_run)
CSV.write("data/stats_per_pair.csv",stats_per_pair)
CSV.write("data/stats_mecamylamine.csv",mecadf)
CSV.write("data/stats_picrotoxin.csv",picrodf)

## Functions for export to javascript
function writeJS(name,variable_name,someDict)
    out_data = JSON.json(someDict)
    open(name, "w") do f
        write(f, "const ",variable_name,"=",out_data)
    end
end

mkpath("js")

## Functions to export fluorescence data (we need dicts for js, and we're selecting the first 6 runs (no drug) for now)
function get_dataDict_per_key(pk,data_dict)
    dat = data_dict[pk][1:min(6,end)]
    if length(dat)==0
        return(0)
    else
        dat = Dict(x[2]["pulseNumber"] => transformAxisArray(x[1]) for x in dat)
    end
end

function transformAxisArray{T}(aa::AxisArray{T,1})
    Dict("x"=>axes(aa,1)[:]-axes(aa,1)[findfirst(axes(aa,1)[:].>3.0)],"y"=>aa)
end

transformAxisArray{T}(aa::AxisArray{T,3}) = Dict("x"=>axes(aa,1)[:]-axes(aa,1)[findfirst(axes(aa,1)[:].>3.0)],
                                                 "y"=>reshape(aa,(size(aa,1)*size(aa,2),size(aa,3))))

transformAxisArray{T}(AA::Array{AxisArrays.AxisArray{T,1,Array{T,1},Tuple{AxisArrays.Axis{:time,StepRangeLen{Float64}}}},1}) =  Dict("x"=>[axes(aa,1)[:]-axes(aa,1)[findfirst(axes(aa,1)[:].>3.0)] for aa in AA],"y"=>AA)

## Exporting the full fluorescence and the average
out_data = Dict(pk => get_dataDict_per_key(pk,full_data_dict) for pk in keyEntries)
filter!((x,y) -> y!=0,out_data)
    
writeJS("js/full_data.js","FULL_DATA",out_data)

out_data_avg = Dict(pk => get_dataDict_per_key(pk,avg_data_dict) for pk in keyEntries)
filter!((x,y) -> y!=0,out_data_avg)
    
writeJS("js/avg_data.js","AVG_DATA",out_data_avg)

## Exporting a table mapping the experiments to the cell pairs
pairToExp = Dict(cpair => convert(Array{String},labbook[convert(Array{Bool},labbook[:cellToCell].==cpair),:keyEntry]) for 
                 cpair in unique(labbook[:cellToCell]))

writeJS("js/pairsToExp.js","PAIRS_TO_EXP",pairToExp)

supertypes = unique(linesToType[:Supertype])

## Exporting a supertype description table
##Limiting to the experiments that have been done
superDict = Dict(s => unique(linesToType[(linesToType[:Supertype].==s) .& 
        [td in uniqueTypesUsed for td in linesToType[Symbol("Type Description")]],Symbol("Type Description")]) for 
    s in supertypes)
superDict = filter((k,x) -> !isempty(x),superDict)

writeJS("js/supertypes.js","SUPERTYPES",superDict)

summaryDataDict = 
Dict(cp => Dict(nP => Dict(string(n) => 
stats_per_pair[(stats_per_pair[:cellPair].==cp) .& (stats_per_pair[:nPulses_median].==nP),
    n][1] for 
    n in names(stats_per_pair)[[2,3,collect(5:end)...]]) for
    nP in stats_per_pair[stats_per_pair[:cellPair].==cp,:nPulses_median]) for 
    cp in unique(stats_per_pair[:cellPair]))
            ## Only exporting the ones that are in the full set (excludes low ns)



## PAIRS REQUIRING OUR ATTENTION FOR REPEATS
superSummary = 
Dict(cp => Dict(string(n) => 
stats_per_pair_20[(stats_per_pair_20[:cellPair].==cp),n][1] for 
            n in [:n,:expType,:signif1,:signif5,:distanceNorm]) for
            cp in unique(stats_per_pair_20[:cellPair]))

## For now we're exporting the stats for the drug free runs
perRunDataDict = 
Dict(exp => Dict(nP => Dict(string(n) => 
unique(stats_per_run[(stats_per_run[:experiment].==exp) .& (stats_per_run[:nPulses_median].==nP) .& (stats_per_run[:preDrug]),
    n]) for 
                n in names(stats_per_run)) for
    nP in stats_per_run[stats_per_run[:experiment].==exp,:nPulses_median]) for 
    exp in unique(stats_per_run[:experiment]))

## A table of drivers, used by the website
drivers = Dict(td => convert(Array,
        linesToType[linesToType[Symbol("Type Description")].== td,:Line]) for 
                          td in unique(linesToType[Symbol("Type Description")]));
writeJS("js/drivers.js","DRIVERS",drivers)

writeJS("js/perRunData.js","PER_RUN_DATA",perRunDataDict)
writeJS("js/neurontypes.js","NEURON_TYPES",neuronTypes)
writeJS("js/summaryData.js","SUMMARY_DATA",summaryDataDict)
writeJS("js/superSummary.js","SUPER_SUMMARY",superSummary)


