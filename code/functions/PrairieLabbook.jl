
### Import the labbook table selecting certain values for a sub-table (returns two dataframes, the full one and the subset).
function readLabbook(tablePath,lines;expDay=nothing,tagExcl=":Problematic:",dateFmt="yyyy-mm-dd e")
    labTable = CSV.read(tablePath,null="NA",dateformat=dateFmt,weakrefstrings=false,quotechar=''')#readtable(tablePath,nastrings=["NA"])

    labTable[:folderName] = lowercase.(Dates.format(labTable[:TIMESTAMP],"uddyy"))
    labTable[:Runs] = collect(vcat(eval(parse(x))...) for x in Missings.replace(labTable[:Runs],"[0]"))
    labTable[:Region] = collect(eval(parse(x)) for x in Missings.replace(labTable[:Region],"Dict(\"Missing\" => [0])"))
    
    labTableFinal = similar(labTable,0)[:,names(labTable) .!= :Region]
    labTableFinal[:Region]="Missing"
    labTableFinal[:RegionRuns] = Missings.missing
    
    for fly in 1:size(labTable)[1]
        for reg in 1:length(labTable[:Region][fly])
            labLine = labTable[fly,names(labTable) .!= :Region]
            labLine[:Region] = collect(keys(labTable[:Region][fly]))[reg]
            labLine[:RegionRuns] = [collect(values(labTable[:Region][fly]))[reg]]
            labTableFinal = vcat(labLine,labTableFinal)
        end
    end

    labTable = labTableFinal

    labTable = labbook_extend(labTable,lines)
    
    if typeof(expDay)==Void
        subsetTable = labTable[find(Missings.replace(labTable[:TAGS].!=tagExcl,true)),:]
    else
        subsetTable = labTable[[in(x,expDay)::Bool for x in labTable[:folderName]] .& collect(Missings.replace(labTable[:TAGS].!=tagExcl,true)),:]
    end
    return((labTable,subsetTable))
end


### Creates a dictionnary containing all the info for a given fly from the corresponding line in the labbook data frame (needs the path to the data because it uses the Prairie xml files)
function makeflyDict(tableLine;dataFolder="/groups/jayaraman/jayaramanlab/DATA/Romain/LALConnectivityProject")

    ### Build the xml file names for the runs covered by this fly.
    dataPathName = "$dataFolder/$(tableLine[:folderName][1])"
    #runs = [eval(parse(tableLine[:Runs][1]))...]

    runStrings = [@sprintf("%03d",r) for r = tableLine[:RegionRuns][1]]
    ## Folder names
    #runTNames = ["/TSeries-$(tableLine[:TIMESTAMP][1])-$(@sprintf("%03d",r))" for r in tableLine[:Runs][1]]
    runTNames = filter(x -> (contains(x,"TSeries") & any([ismatch(Regex("$y\$"),x) for y=runStrings])),readdir(dataPathName))

    ## Full path names
    runNames = pmap(runTNames) do runN
        fileNs = readdir("$dataPathName/$runN")
        fileN = filter(x -> (contains(x,".xml") & !contains(x,"Voltage")), fileNs)[1]
        "$dataPathName/$runN/$fileN"
    end

    map(importPrairie,runNames)
end

### Computes the index of the fly in a given experimental day
function getFlyN(tableLine,labTable)
    findfirst(labTable[labTable[:folderName].==tableLine[:folderName],:][:Runs].==tableLine[:Runs])
end



### Processing for a single run

function runExtract(prairie;mvtCorrect = false,repeatsReg = false,channel=2)
    nEp = length(prairie["sequences"])                    ## Number of repeats in the run

    ### Import the data per se
    greenCh = map(1:nEp) do i
        getPrairieFrames(prairie,channel=channel,seqN=i)#*1.0
    end

    ### Movement correction of the individual repeats
    if mvtCorrect
        for (i,im) in enumerate(greenCh)
            greenCh[i] = alignFromDict(im,stackDftReg(im))
        end
    end

    ### Averages of each repeat
    #greenChAvs = pmap((x) -> getindexim(meanfinite(x,3),:,:,1),greenCh)
    greenChAvs = pmap((x) -> squeeze(Images.meanfinite(x.data,AxisArrays.axisdim(x,AxisArrays.Axis{:time})),AxisArrays.axisdim(x,AxisArrays.Axis{:time})),greenCh)

    ### Global movement correction of the repeats (aligned on the first repeat)
    if repeatsReg
        repeatsDft = Array(Array{Dict},length(greenChAvs)-1)

        for i in 2:length(greenChAvs)
            repeatsDft[i-1] = stackDftReg(greenChAvs[i],ref=greenChAvs[1])
            greenChAvs[i] = alignFromDft(greenChAvs[i],repeatsDft[i-1])
            toShift = repeatsDft[i-1]
            toShift["shift"] = [toShift["shift"];0] ## This is to do the registration on 3d images, with no shift in z
            greenCh[i] = alignFromDft(greenCh[i],toShift)
        end
    end

    ## Rerun the indiviual movement corrections with the average image as a reference
    if mvtCorrect
        for (i,im) in enumerate(greenCh)
            greenCh[i] = alignFromDft(im,stackDftReg(im,ref=greenChAvs[i]))
        end
        greenChAvs = map((x) -> meanfinite(x.data,axisdim(x,Axis{:time})),greenCh)
    end

    ## Full average
    greenChAv = mean(greenChAvs)

    ### Get a rough estimate of the background from that image.
    max4B = round(Integer,length(greenChAv)/10)
    B = mean(sort(vec(greenChAv))[1:max4B])

    ## Average response
    runAv = mean(greenCh)

    ## DAC stim, considering only pulse trains for now
    #assuming same pulse train for all repeats
    if (prairie["DAC"][1] == nothing)
        pulseTrain = [3]
        stimIntensity = 1
    else
        pulseTrain = dac2Train(prairie["DAC"][1][1])
        if typeof(prairie["DAC"][1][1]) == String
            dac = prairie["DAC"][1][1]
            prm = open(dac)
            dac = readdlm(prm,'=','\n')
            stimIntensity = dac[dac[:,1].== "DAC1 Train Pulse Pot 1",2]/1000
        else
            stimIntensity = prairie["DAC"][1][1]["DACStimulusParameters"]["PulsePotentialStart"]
        end
    end
    pulseNumber = length(pulseTrain)
    runStart = DateTime(prairie["sequences"][1]["attributes"]["time"][1:8],"H:M:S")
    runStop = DateTime(prairie["sequences"][end]["attributes"]["time"][1:8],"H:M:S")+Dates.Second(20)


    Dict("green" => greenCh,"av"=> greenChAv, "B" => B, "runAv" => runAv, "pulseNumber"=>pulseNumber,"stimIntensity"=>stimIntensity,"runStart"=>runStart,"runStop"=>runStop)
end



### Very primitive, only works for a single pulse train on DAC1
function dac2Train(dac)
    if typeof(dac) != String
        startPt = dac["DACDeviceParameters"]["relativeTime"]+dac["DACStimulusParameters"]["FirstPulseDelay"]
        inter = dac["DACStimulusParameters"]["PulseWidth"]+dac["DACStimulusParameters"]["PulseSpacing"]
        nPul = dac["DACStimulusParameters"]["PulseCount"]
    else
        prm = open(dac)
        dac = readdlm(prm,'=','\n')
        startPt = dac[dac[:,1].== "DAC1 Time Delay to First Pulse Train",2][1] + dac[dac[:,1].=="DAC1 Train First Pulse Delay 1",2][1]
        inter = dac[dac[:,1].== "DAC1 Train Pulse Duration 1",2][1] + dac[dac[:,1].== "DAC1 Train Inter-Pulse Delay 1",2][1]
        nPul = int(dac[dac[:,1].== "DAC1 Train Num Pulses 1",2][1])
    end
    collect(range(startPt/1000,inter/1000,nPul))
end


### Add some columns to the labbook using a table describing the lines and types
function labbook_extend(labbook_table,linesToType)
    sort!(labbook_table,cols=[:folderName])
    labbook_table[:genotype] = labbook_table[:LexA].*"LexA-".*labbook_table[:Gal4].*
    "Gal4-Chrimson-in-".*labbook_table[:ActivatorExp]
    labbook_table[:genotypeRegion] = labbook_table[:genotype] .*"-".* labbook_table[:Region]
    labbook_table[:genotypePre] = [labbook_table[Symbol(labbook_table[:ActivatorExp][i])][i]*
                                   labbook_table[:ActivatorExp][i] for 
                                   i in 1:size(labbook_table,1)]
    labbook_table[:genotypePost] =[(labbook_table[:ActivatorExp][i] == "Gal4" ? 
                                    labbook_table[:LexA][i] : labbook_table[:Gal4][i])*
                                   labbook_table[:ReporterExp][i] for i in 1:size(labbook_table,1)]
    labbook_table[:driverPre] = [labbook_table[Symbol(labbook_table[:ActivatorExp][i])][i] for 
                                 i in 1:size(labbook_table,1)]
    labbook_table[:driverPost] =[labbook_table[:ActivatorExp][i] == "Gal4" ? 
                                 labbook_table[:LexA][i] : labbook_table[:Gal4][i] for 
                                 i in 1:size(labbook_table,1)]
    labbook_table[:cellPre] = [linesToType[Symbol("Type Description")][
                                                              linesToType[:Line] .== labbook_table[:driverPre][i]][1]  for
                                                          i in 1:size(labbook_table,1)]
    labbook_table[:cellPost] = [linesToType[Symbol("Type Description")][
                                                               linesToType[:Line] .== labbook_table[:driverPost][i]][1] 
                                for i in 1:size(labbook_table,1)]
    labbook_table[:cellToCell] = labbook_table[:cellPre] .*"-to-".* labbook_table[:cellPost]
    ## We use the unique run numbers to identidy the flies, as several regions might have 
    ## been imaged for the same fly.
    flyN = by(labbook_table,:folderName) do df
        DataFrame(flyN=[findfirst(df[:Runs],ru) for ru in df[:Runs]])
    end
    
    labbook_table[:flyN]=flyN[:flyN]
    labbook_table[:keyEntry] = labbook_table[:folderName].*"-Fly".*[string(labbook_table[i,:flyN]) for i in 1:size(labbook_table)[1]].*"-".*labbook_table[:Region]
    
    
    labbook_table
end
