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
    sort!(mini_labbook,:keyEntry)
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
    sort!(mini_labbook,:keyEntry)
    
    pair_keys = mini_stats[:experiment] 
    mini_data = [get(data_dict,pair_keys[i],0)[mini_stats[:runIdx][i]][1] for 
                                                         i in eachindex(pair_keys)];
    Dict("stats"=>mini_stats,"data"=>mini_data)
end

function select_data(cellPair::String,nPulses,data_dict,run_range,labbook,stats_per_run=stats_per_run)
    select_data([cellPair],nPulses,data_dict,run_range,labbook,stats_per_run)
end

### Plot using the output of the select_data function
function make_raw_plot(dataset;substract=false,scalebar=true,scaley=3,np=20,colorV=:experiment,legendV=false,traceW=3,kwargs...) ## substract is a keyword argument deciding if the
                                                          ## baseline fluorescence is to be substracted for display
    p = plot(textcolor=RGB(128/255,128/255,128/255),guidecolor=RGB(128/255,128/255,128/255);kwargs...)
    make_raw_plot!(p,dataset;substract=substract,scalebar=scalebar,scaley=scaley,np=np,colorV=colorV,legendV=legendV,traceW=traceW,kwargs...)
    p
end

function make_raw_plot!(p,dataset;substract=false,scalebar=true,scaley=3,np=20,colorV=:experiment,legendV=false,pairTitle=false,traceW=3,kwargs...) ## substract is a keyword argument deciding if the
                                                          ## baseline fluorescence is to be substracted for display
    topvalue = 0
    colorsVIdx = unique(dataset["stats"][colorV])

    if pairTitle
        plot!(p,title=dataset["stats"][:shortPair][1])
    end
 
        
    for i in 1:length(dataset["data"])
        d = copy(dataset["data"][i])
        if substract 
            d[:] = d .- dataset["stats"][:baseline_median][i]
        end
        plot!(p,d,color=find(colorsVIdx.==dataset["stats"][colorV][i])[1],alpha=0.8,
            lw=traceW,label="";kwargs...)
        topvalue = max(maximum(d),topvalue)
    end

    
    if legendV
        colorsVI = [findfirst(dataset["stats"][colorV],g) for g in colorsVIdx]
        println(colorsVIdx)
        for i in 1:length(colorsVI)
            p.subplots[1].series_list[colorsVI[i]][:label] = colorsVIdx[i]
        end
    end
    plot!(p,[0;0.033*np],[topvalue;topvalue],color=:gray80,fill=:gray80,alpha=0.4,fillrange=0,line=:path,label="";kwargs...)
    if scalebar
        plot!(p,[6;8],[scaley-1;scaley-1],color=:gray50,line=:path,lw=3,label="";kwargs...)
    end
end

## Histograms segregated by experiment type
function makeStatHist(statDf,statN;kwargs...)
        stephist(statDf[statN],group=statDf[:expType],fill=true,alpha=0.4,normalize=:none,
            nbins=linspace(minimum(statDf[statN]),maximum(statDf[statN]),30),
                 title=fullNamesDF[statN][1],titlefontsize=10
                 ;kwargs...)
end

function makeStatHist!(statDf,pl,stat;kwargs...)
    stephist!(pl,statDf[stat],group=statDf[:expType],fill=true,alpha=0.4,normalize=:none,
            nbins=linspace(minimum(statDf[stat]),maximum(statDf[stat]),30),
              xlab=""
            ;kwargs...)
end

## Matrix plot

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

function makeMatrixPlot(statUsed,namesUsed;kwargs...)
    mat = statsMatrices[statUsed]
    matVals = filter(x->!isnan(x),mat)
    grad = ColorGradient(ColorGradient(:bluesreds).colors,[0,-minimum(matVals)/(maximum(matVals)-minimum(matVals)),1.0])

    themat = heatmap(namesUsed, namesUsed, mat,aspect_ratio=1,
                     color=grad,
                     xrotation=90
                     ,yguide="Presynaptic candidate",xguide="Postsynaptic candidate",
                     xticks = (0.5:1:length(namesUsed)-0.5,namesUsed),
                     yticks = (0.5:1:length(namesUsed)-0.5,namesUsed);
                     kwargs...
         )
    plot!(themat,identity,line=(:grey80,:dot),lab="")
    scatter!(themat,collect(findn(matGuesses))-0.5...,m=(:black),msw=0,malpha=0.5,lab="Overlapping pairs")
    themat
end

#function makePairDrugPlots(df,cp;kwargs...)
#    pairPlot = plot(layout=(1,2);kwargs...)
#    makePairDrugPlots!(pairPlot,df,cp,1)
#end
# Drug summaries
function makePairDrugPlots(df,cp;drugTime=[6;11])
    preData = []
    drugData = []
    postData = []
    df[:isWash] = false
    for expe in unique(df[df[:cellPair].==cp,:experiment])
        subDf = df[df[:experiment].==expe,:]
        df[find(df[:experiment].==expe)[(size(subDf,1)-1):size(subDf,1)],:isWash]=true
        runPre = subDf[subDf[:timeToDrug].<=2,:runIdx]
        runDrug = subDf[(drugTime[1].<subDf[:timeToDrug].<drugTime[2]),:runIdx]
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

    # In case there's only one experiment
    preDataS[isnan.(preDataS)]=0
    drugDataS[isnan.(drugDataS)]=0
    postDataS[isnan.(postDataS)]=0
    

    topvalue = maximum([maximum(preDataM.+preDataS),maximum(postDataM.+postDataS),maximum(drugDataM.+drugDataS)])+0.2
    

    ystart = max(minimum(df[df[:cellPair].==cp,:integNorm_scaled]),-1.0)-0.05
    yend = min(maximum(df[df[:cellPair].==cp,:integNorm_scaled]),1.0)+0.05

    plotTimeCourse = @> df[df[:cellPair].==cp,:] begin
                        @splitby (_.experiment,_.genotype,_.shortPair)
                        @x _.timeToDrug
                        @y _.integNorm_scaled#_.integral_to_peak_scaled #_.distanceNorm 
                        @set_attr :title _[3]
                        @set_attr :hover _[2]
                        # @set_attr :color _[4] ? 3 : :gray40
                        @plot plot(legend=:none,msw=0,label="",color=:gray40,title_location=:right,width=3,alpha=0.6,link=:x,ylim=(ystart,yend))
                        end
    
 @> df[(df[:cellPair].==cp) .& df[:isWash],:] begin
        @splitby (_.experiment,_.genotype,_.shortPair,_.isWash)
        @x _.timeToDrug
        @y _.integNorm_scaled#_.integral_to_peak_scaled #_.distanceNorm 
        @plot plot!(plotTimeCourse,legend=:none,
                    msw=0,label="",color=3,title_location=:right,width=3,alpha=0.6,link=:x,ylim=(ystart,yend))
    end

    plotAv = plot(preDataM,ribbon=preDataS,label="",lw=3,top_margin=7mm,bottom_margin=5mm,yaxis=false,link=:x)
    plot!(plotAv,drugDataM,ribbon=drugDataS,label="",lw=3)
    plot!(plotAv,postDataM,ribbon=postDataS,label="",lw=3)

    plot!(plotAv,[0;0.033*np],[topvalue;topvalue],color=:gray80,fill=:gray80,alpha=0.4,
        fillrange=0,line=:path,label="")

    plot!(plotAv,[-2;8],[0;0],color=:gray80,line=(3,:dash),label="")

    plot!(plotTimeCourse,[0;0],[ystart;yend],color=:gray40,line=(3,:dash),label="")

    ## Adding shading for the different periods
    ## Control
    bar!(plotTimeCourse,[-4;0],[yend;yend],color=1,fill=1,alpha=0.3,
        line=:path,label="")
    plot!(plotTimeCourse,[-4;0],[ystart;ystart],color=1,fill=1,alpha=0.3,
        line=:path,label="")

    ## Drug
    plot!(plotTimeCourse,[drugTime[1];drugTime[2]],[ystart;ystart],color=2,fill=2,alpha=0.3,
        line=:path,label="")
    plot!(plotTimeCourse,[drugTime[1];drugTime[2]],[yend;yend],color=2,fill=2,alpha=0.3,
        line=:path,label="")

    ## Scale Bars
    plot!(plotAv,[9.3,9.3],[0;1],color=:gray40,lw=4,label="")
    
    
    [plotTimeCourse,plotAv]
end
