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

### Plot using the output of the select_data function
function make_raw_plot(dataset;substract=false,scalebar=true,scaley=3,np=20,colorV=:experiment,traceW=3,kwargs...) ## substract is a keyword argument deciding if the
                                                          ## baseline fluorescence is to be substracted for display
    p = plot(textcolor=RGB(128/255,128/255,128/255),guidecolor=RGB(128/255,128/255,128/255);kwargs...)
    make_raw_plot!(p,dataset;substract=substract,scalebar=scalebar,scaley=scaley,np=np,colorV=colorV,traceW=traceW,kwargs...)
    p
end

function make_raw_plot!(p,dataset;substract=false,scalebar=true,scaley=3,np=20,colorV=:experiment,traceW=3,kwargs...) ## substract is a keyword argument deciding if the
                                                          ## baseline fluorescence is to be substracted for display
    topvalue = 0
    colorsVIdx = unique(dataset["stats"][colorV])
    for i in 1:length(dataset["data"])
        d = copy(dataset["data"][i])
        if substract 
            d[:] = d .- dataset["stats"][:baseline_median][i]
        end
        plot!(p,d,color=find(colorsVIdx.==dataset["stats"][colorV][i])[1],legend=false,alpha=0.8,
            lw=traceW;kwargs...)
        topvalue = max(maximum(d),topvalue)
    end
    plot!(p,[0;0.033*np],[topvalue+1;topvalue+1],color=:gray80,fill=:gray80,alpha=0.4,fillrange=0,line=:path;kwargs...)
    if scalebar
        plot!(p,[6;8],[scaley-1;scaley-1],color=:gray50,line=:path,lw=3;kwargs...)
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

function makeMatrixPlot(statUsed;kwargs...)
    mat = statsMatrices[statUsed]
    matVals = filter(x->!isnan(x),mat)
    grad = ColorGradient(ColorGradient(:bluesreds).colors,[0,-minimum(matVals)/(maximum(matVals)-minimum(matVals)),1.0])

    themat = heatmap(uniqueTypesUsed, uniqueTypesUsed, mat,aspect_ratio=1,
                     color=grad,
                     xrotation=90
                     ,yguide="Presynaptic candidate",xguide="Postsynaptic candidate",
                     xticks = (0.5:1:length(uniqueTypesUsed)-0.5,uniqueTypesUsed),
                     yticks = (0.5:1:length(uniqueTypesUsed)-0.5,uniqueTypesUsed);
                     kwargs...
         )
    plot!(themat,identity,line=(:grey80,:dot),lab="")
    scatter!(themat,collect(findn(matGuesses))-0.5...,m=(:black),msw=0,malpha=0.5,lab="Overlapping pairs")
    themat
end

function makePairDrugPlots(df,cp)
    pairPlot = plot(layout=(1,2))
    makePairDrugPlots!(pairPlot,df,cp,1)
end
# Drug summaries
function makePairDrugPlots!(pairPlot,df,cp,subStart)
    preData = []
    drugData = []
    postData = []
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

    # In case there's only one experiment
    preDataS[isnan.(preDataS)]=0
    drugDataS[isnan.(drugDataS)]=0
    postDataS[isnan.(postDataS)]=0
    
    plot!(preDataM,ribbon=preDataS,subplot=subStart+1,xlabel="Time (s)",label="",right_margin=10mm,top_margin=10mm,bottom_margin=10mm,lw=3,yticks=4)
    plot!(drugDataM,ribbon=drugDataS,subplot=subStart+1,label="",lw=3,yticks=4)
    plot!(postDataM,ribbon=postDataS,subplot=subStart+1,label="",lw=3,yticks=4)
    
    topvalue = maximum([maximum(preDataM.+preDataS),maximum(postDataM.+postDataS),maximum(drugDataM.+drugDataS)])+1
    plot!([0;0.033*np],[topvalue;topvalue],color=:gray80,fill=:gray80,alpha=0.4,
        fillrange=0,line=:path,subplot=1+subStart,label="")

    ystart = max(minimum(df[df[:cellPair].==cp,:integNorm_scaled]),-1.0)-0.1
    yend = min(maximum(df[df[:cellPair].==cp,:integNorm_scaled]),1.0)+0.1
    
    @> df[df[:cellPair].==cp,:] begin
        @splitby (_.experiment,_.genotype,_.shortPair)
        @x _.timeToDrug
        @y _.integNorm_scaled#_.integral_to_peak_scaled #_.distanceNorm 
        @set_attr :title _[3]
        @set_attr :hover _[2]
        @plot scatter!(legend=:none,
                       ylabel="Response integral",xlabel="Time to drug",subplot=subStart,msw=0,ylim=(ystart,yend),label="",titlefont=font("DejaVu Sans",8),title_location=:right,yticks=4)
    end
    pairPlot
end
