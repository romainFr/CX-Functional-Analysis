
## This script reads the data, detects the ROIs and returns the results in a
## Dictionnary. To be used on the command line (ideally in parallel) with "julia -p X functionalConnectivityRaw.jl path_to_data_folder nameOfFolder1ToAnalyze nameOfFolder2ToAnalyze ..."
## If arguments 2 and up are omitted, all experiments are analyzed.

using JLD2
using DataFrames,CSV
using Images
using FileIO,ImageMagick
using SubpixelRegistration
using FluorescentSeries
using PrairieIO
using Clustering

include("functions/PrairieLabbook.jl")
include("functions/imageUtilities.jl")

############################################## Parameters definition ###########
## Where the labbook table is :
tablePath = "labbookTable.csv"

### If true, just add the processed entry to the existing table, otherwise creates a new table from scratch
toUpdate = false

### Where the data is :
baseDataFolder =  ARGS[1]#"../dmHere/DATA/Romain/LALConnectivityProject/"
                         #"example-raw-data/"
### Which experimental day we want to process (or nothing if we want to process everything), this is usually passed by argument.
#expday=nothing ["feb1916"]
if length(ARGS) == 1
    expday = nothing
else
    expday = ARGS[2:end]
end
println(expday)

try
    mkdir("data")
end
#################################################################################

problemFolders = String[]

linesToType = CSV.read("LinesAndTypes.csv")
(mainTab,subTab) = readLabbook(tablePath,linesToType,expDay=expday)

if toUpdate == false
    fullDf = Dict()
    images = Dict()
else
    fullDf = load("data/rawData.jld2")
    images = load("data/expImages.jld2")
end

for i in (1:size(subTab)[1])
    fly = makeflyDict(subTab[i,:],dataFolder=baseDataFolder);

    pockelsParams = [fly[i]["globalConfig"]["laserPower_0"] for i in eachindex(fly)]

    if (!ismissing(subTab[i,:][:Drug][1]))
        drugStart = DateTime(subTab[i,:][:DrugTime][1][1:8],"H:M:S")
    else
        drugStart = Missings.missing
    end

    info(subTab[i,:keyEntry])

    if length(fly) == 0
        println("Nothing in there")
        push!(problemFolders,subTab[i,:keyEntry])
        continue
    end
    ## Load all the data to get the baseline and the ROIs
    info("Load data")
    runFluos = map(fly) do run
        gc()
        roiDict = runExtract(run,mvtCorrect=false)
    end
    info("Loaded")
    
    ## Checking for inhomogenous ROIs
    if any(runFluos .== "Run not completed")
        println("Bad runs selection, one couldn't be loaded")
        push!(problemFolders,subTab[i,:keyEntry])
        continue
    end

    if any([size(runFluos[run]["av"]) != (size(runFluos[1]["av"])) for run in 2:length(runFluos)])
        println("Bad runs selection, some have a different size...")
        push!(problemFolders,subTab[i,:keyEntry])
        continue
    end
    ## Alignment of the different runs, relatively low resolution (so pixel statistics are not affected)
    info("Align runs")

    shifts = SharedArray{Int64}(2,length(runFluos)-1)
    ref = convert(SharedArray,runFluos[length(runFluos)]["av"])

    runFluos[1:(length(runFluos)-1)] = pmap(runFluos[1:(length(runFluos)-1)],1:(length(runFluos)-1)) do rF,run
        registration = subpixelRegistration.stackDftReg(rF["av"],ref=ref,ufac=1)
        shifts[:,run] = round.(Int64,registration["shift"])
        ## Align the grand averages
        rF["av"] = subpixelRegistration.alignFromDict(rF["av"],registration)
        ## Align the run averages (those are "volumes" so we need to set the z shift to 0)
        registration["shift"] = [registration["shift"];0]
        rF["runAv"] = subpixelRegistration.alignFromDict(rF["runAv"],registration)
        ## Finally align the individual runs
        for rep in 1:length(rF["green"])
            rF["green"][rep][:,:,:] = subpixelRegistration.alignFromDict(rF["green"][rep],registration)
        end
        rF
    end
    maxShiftPos = round.(Int64,maximum([[0;0] shifts],2))
    maxShiftNeg = round.(Int64,minimum([[0;0] shifts],2))
    
    ## Cropping everything to the area common to all repeats // catch time to drug
    cropRegion = [(1+maxShiftPos[i]):(size(runFluos[1]["av"])[i]+maxShiftNeg[i]) for i in 1:2]
    for run in 1:length(runFluos)
        runFluos[run]["av"] = runFluos[run]["av"][cropRegion...]
        runFluos[run]["runAv"] = runFluos[run]["runAv"][cropRegion...,:]
        for rep in 1:length(runFluos[run]["green"])
            runFluos[run]["green"][rep] = runFluos[run]["green"][rep][cropRegion...,:]
        end
        if !ismissing(drugStart)
                runFluos[run]["timeToDrug"] = drugStart - runFluos[run]["runStart"]
        else
            runFluos[run]["timeToDrug"] = Missings.missing
        end
    end
    info("Aligned")
    ## Average of all the repeats
    grdAv = mapreduce(x -> x["av"],+,runFluos)/length(runFluos)

    info("Get ROIs")

    fluorescentRegions = kMeansIm(grdAv,2)
    
    info("Compute fluorescence")
    fluoSimple = [cat(3,[FluorescentSerie(rep,fluorescentRegions) for rep in runF["green"]]...) for runF in runFluos]
    fluoParams = [Dict("pulseNumber" => runF["pulseNumber"],"timeToDrug" => runF["timeToDrug"],"stimIntensity"=>runF["stimIntensity"]) for runF in runFluos]
    ## Find the global F0 as the median of the 3% dimmer pixels ?

    deltaFluoSimple = fluoSimple
    globalBackground = mean([run["B"] for run in runFluos])
    powerLevels = unique(pockelsParams)
    if length(powerLevels)>1
      info("More than one Pockels setting used here")
    end

    globalBaseline = Array{Array{Float64,1}}(length(powerLevels))
    for pw in eachindex(powerLevels)
        runsIdx = find(pockelsParams .== powerLevels[pw])
        if length(runsIdx) == 1
            giantFluo = fluoSimple[runsIdx[1]]
        else
            giantFluo = vcat([x.data for x in fluoSimple[runsIdx]]...)
        end
        giantFluo = reshape(permutedims(giantFluo,[1;3;2]),size(giantFluo,1)*size(giantFluo,3),size(giantFluo,2))
        globalQ3 = [quantile(giantFluo[:,i],0.05) for i in 1:size(giantFluo)[2]]
        globalBaseline = median(giantFluo[find(giantFluo.< globalQ3.')],1)[:]
        deltaFluoSimple[runsIdx] = [deltaFF(fluoSimple[runsIdx][i],globalBaseline,globalBackground) for i in eachindex(fluoSimple[runsIdx])]
    end

    fullDf[subTab[i,:keyEntry]] = collect(zip(deltaFluoSimple,fluoParams))

    grdAv = AxisArray(Gray.(N4f12.(grdAv)),axes(runFluos[1]["green"][1],1),axes(runFluos[1]["green"][1],2))
    fluorescentRegions = AxisArray(fluorescentRegions,axes(runFluos[1]["green"][1],1),axes(runFluos[1]["green"][1],2))

    images[subTab[i,:keyEntry]] = Dict("average_image"=> grdAv,"ROIs"=> fluorescentRegions)
end

## Adding some columns to the labbook and cleaning it for further use
info("Export results")
save("data/rawData.jld2",fullDf)
save("data/expImages.jld2",images)

info("Export labbook")
mainTab = mainTab[ismissing.(mainTab[:TAGS]),:]
mainTab[:timesToDrug]=Array{Array,1}(size(mainTab,1))

for k in keys(fullDf)
     mainTab[findfirst(mainTab[:keyEntry].==k),:timesToDrug] = [kdct[2]["timeToDrug"] for kdct in fullDf[k]]
end

save("data/labbookTable.jld2","labbook",mainTab)
