## A very simple k-means applied to image intensity
function kMeansIm(img,nl)
    imgR = reshape(img,(size(img)[1]*size(img)[2],nimages(img))).'
    kRes = kmeans(imgR,nl)
    ### We want the cluster numbers to be sorted by the average intensity in the cluster
    reord = sortperm(mean(kRes.centers,1)[:])
    clusts = map((x) -> findfirst(reord,x),assignments(kRes))
    roi = reshape(clusts,size_spatial(img))-1
    roi
end
