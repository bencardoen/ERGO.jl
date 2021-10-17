# This file is part of ERGO.jl.
#
#  ERGO.jl is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License v3 as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ERGO.jl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License v3 for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ERGO.jl.  If not, see <https://www.gnu.org/licenses/>.
# Copyright Ben Cardoen, 2018-2021

module ERGO

import DataFrames
import CSV
using Logging
import Images
import Random
import Glob
import Statistics
import JSON
using LinearAlgebra
import ProgressMeter.@showprogress


export parseFrame, writeTrainingData, gmsm, gm, check, to_cc, harmonic_pair, geom_pair, peaks, centernorm,
buildAdjMat, scoreCC, framecorrelation, parseFixedImages, parseFixedFrame, writeTrainingData,
findCentroids, buildGraph, buildSparseGraph, connectedComponents, collect_results, procresults,
centernorm, normimg, tomask, aszero, mcc, pearsoncorrelation


"""
    tomask(img)

Utility function to binarize argument
"""
function tomask(img)
    c = copy(img)
    c[c .> zero(eltype(img))] .= oneunit(eltype(img))
    return Images.Gray{Images.N0f8}.(c)
end

"""
	normimg(img[, center=false])

Normalize 2/3D image
if center = true : (a+min(a))/range(a)
"""
function normimg(img, center=false)
    if iszero(img)
		return copy(img)
	end
	@assert(! any(isnan.(img)))
	if length(size(img)) == 2
		return _normimg(img, center)
	end
	if length(size(img)) == 3
		res = aszero(img)
		for z in 1:size(img, 3)
			res[:,:,z] = _normimg(img[:,:,z],center)
		end
		return res
	end
    @error "Normalizing image with dims not in [2,3]"
	@assert false
end

function _normimg(img, center=false)
	if center == true
		return centernorm(img)
	end
	return img ./ maximum(img)
end

"""
	aszero(x)

Utility function to return a zeroed copy of the argument
"""
function aszero(x)
    return zeros(eltype(x), size(x))
end

"""
	Report results in a dataframe
	outdir : path to save, file will be saved as outdir/intensity_kn.csv
	kn : postfix to prevent overwriting existing files
	res : array of results, length(array) >= fcount
	fcount : only save results for 1:fcount frames
	Returns an array of fcountx1 values, where array[i] is the number of local maxima in frame i with intensity exceeding μ + kσ (see function parseFixedFrame)
"""
function procresults(res, outdir, kn, fcount)
    N = 0
    M = 0
    pks = zeros(Float64, fcount)
    intensity_data = DataFrame(framenr=Int[], mean=Float64[], std=Float64[],threshold=Float64[],score=Float64[],distance=Float64[])
    for i in 1:fcount
        img, d, score, CC, weightedcentroid, ROIS, modifiedcenter, NegBox, minfiltered, maxfiltered, t = res[i];
        N += size(maxfiltered,1)
        _, _, μ, σ, intens = peaks(img)
        push!(intensity_data, [i, μ, σ, t, score, d])
        pks[i] = size(maxfiltered,1)
    end
    CSV.write(joinpath(outdir,"intensity_$kn.csv"), intensity_data);
    return pks
end

"""
	Writes a csv (outdir/density_label.csv) that can be used for training.
	Information for frame i is stored in results[i], tdata is a training data dict.
	The final csv contains following columns
	framenr roinr mx my count
	Where
		- framenr : > 0 frame index
		- roinr : reference the ROI in order of processing, 0 indicates a 'negative' ROI, e.g. where no emission is assumed
		- count : nr of emission in ROI
		- mx, my : lower left corner of ROI in indices (in frame)
"""
function collect_results(tdata, results, outdir, label="")
    outdf = DataFrame(framenr=Int[], roinr=Int[], mx=Int[], my=Int[],count=Int[])
    for (frame_nr, values) in tdata["data"]
        _, _, _, _, _, ROIS, _, NegBox, _, _, t = results[frame_nr]
        for (roi_nr, rvalues) in values
            mx, my = 0, 0
            if roi_nr == "missed"
                continue
            end
            if roi_nr == -1
                roi_nr = 0
                x1, x2 = NegBox[1][1][1], NegBox[1][2][1]
                y1, y2 = NegBox[1][1][2], NegBox[1][2][2]
                mx = min(x1, x2)
                my = min(y1, y2)
            else
                x1, x2 = ROIS[roi_nr][1,1], ROIS[roi_nr][2,1]
                y1, y2 = ROIS[roi_nr][1,2], ROIS[roi_nr][2,2]
                mx = min(x1, x2)
                my = min(y1, y2)
            end
            N = size(rvalues["roi"],1)
            push!(outdf, [frame_nr, roi_nr, mx, my, N])
        end
    end
    CSV.write(joinpath(outdir,"density_$(label).csv"), outdf)
end

"""
    gmsm(vals)

Compute the geometric mean and standard deviation. Only meaningful for values R+
Vals is filtered for zero entries. (vals = vals[vals .> zero(eltype(vals))]) to avoid collapsing results.

Returns geometric {mean, standard deviation}

"""
function gmsm(vs)
	vals = copy(vs)
    vals = vals[vals .> zero(eltype(vals))]
    N = length(vals)
    li = [log(i) for i in vals]
    gm = exp(sum(li) / N)
    gs = exp(sqrt(sum([log(i/gm)^2 for i in vals])/N))
    return gm, gs
end

"""
	gm(vals)

Return the geometric mean of non-zero values
"""
function gm(vs)
	return gmsm(vs)[1]
end
"""
	From a set of 2D indices (integer) with pixel dimension px, compute connected
	components based on a graph representation (~ spanning tree)
"""
function to_cc(selected, px)
    @info "Building vertices"
    vertices = [[i[1], i[2]] for i in selected];
    distance = sqrt(px^2+px^2)
    @info "Building sparse graph"
    graph = buildSparseGraph(vertices, distance)
    @info "CC"
    CC = connectedComponents(size(vertices)[1], graph)
    @info "Filtering"
    filtered_cc = [cc for cc in CC if length(cc)>2 ];
    return filtered_cc
end

"""
	Harmonic mean for tuple (of arrays)
"""
function harmonic_pair(x, y)
    res = geom_pair(x,y)
    xs = x[res.!=0]
    ys = y[res.!=0]
    res[res.!=0] = (2 .* xs .* ys) ./ (xs .+ ys)
    return res
end

"""
	Geometric mean for tuple
"""
function geom_pair(x, y)
    return sqrt.(x.*y)
end

"""
	Shorthand to generate an array of images from an array of image file names
"""
function loadimages(imagenames)
    return _loadimages(imagenames) |> collect
end

"""
	Shorthand to create a generator of images from an array of image file names.
	E.g. does not load until you iterate over it.
"""
function _loadimages(imagenames)
    return (Images.load(imagename) for imagename in imagenames)
end

"""
	Frame correlation (Pearson)
"""
function pearsoncorrelation(imgi, imgj)
    μ_i = Statistics.mean(imgi)
    μ_j = Statistics.mean(imgj)
    rij = sum((imgi .- μ_i).*(imgj .- μ_j)) / (√( sum((imgi .- μ_i).^2)) * √( sum((imgj .- μ_j).^2)) )
    return rij
end

"""
	For a sequence of images, compute the step separated frame correlation, e.g. between image i and i+step.
	If randomize, first shuffle the image sequence.
	If randomize == true, and seed is not -1, use seed to set RNG.
	Returns array of 1x|images| of correlation values
"""
function framecorrelation(images, randomize=false, seed=0, step=1)
    N = size(images, 1)
    @assert(0 <= step < N-1)
    if randomize
		if seed != -1
        	Random.seed!(seed)
			@debug "Using seed $seed for RNG"
		else
			@warn "Seed not set, not reproducible results will follow !!"
		end
        images = Random.shuffle(images)
	end
    return [pearsoncorrelation(images[i], images[i+step]) for i in 1:N-step]
end

"""
	Finds local minima, maxima, μ, σ and flattened intensity values of img.
"""
function peaks(img)
    intens = Float32.(img[:])
    μ,σ  = Statistics.mean(intens), Statistics.std(intens)
    lM = Images.findlocalmaxima(img);
    lm = Images.findlocalminima(img);
    return lm, lM, μ, σ, intens
end

"""
	Process a set of images
	path : directory where images are stored, with filename [0-9]*.tif, e.g. 001.tif is loaded, abc.tif is not.
	ROI_PX : dimension of region of interest to be extracted, e.g. ROI_PX = 2 --> (2*2+1)^2 ROI. Can be set in function of PSF FWHM, e.g. should be at least large enough to encompass 1 PSF's FWHM.
	range : find only peaks where intensity > μ + kσ , k in range
	train : if true, expects groundtruth to be a dataframe with known GT
	see function parseFixedFrame
	TODO : threaded (or channels)
"""
function parseFixedImages(path, ROI_PX, PX_NM, groundtruth=nothing, range=(1,3), train=false, fixed=0, step=0.1)
    tiffiles = Glob.glob("[0-9]*.tif",path)
    N = size(tiffiles,1)
    @assert(N > 0)
	@debug "Have $N images to process"
    results = Dict()
    for (fi,filename) in enumerate(tiffiles)
        m = match(r"\d+.tif", filename)
        @assert(m != nothing)
        FRAME_NR = parse(Int, filename[match(r"\d+.tif", filename).offset:(end-4)])
        @assert(FRAME_NR > 0)
        results[FRAME_NR] = parseFixedFrame(filename, ROI_PX, PX_NM, FRAME_NR, groundtruth, range, train, fixed)
    end
    return results
end

"""
	Decompose an image (frame) into regions of interest centered on suspected emitters.
	path: file path
	FRAME_NR: used to look up gt (if given)
	groundtruth_dataframe: dataframe with columns frame, xnano, ynano describing emitter locations
	range: parameter sweep of drop  (intensity < μ + k * σ)
	fixed: set intensity threshold manually (if not zero, else find it per frame, in μ+kσ < X < μ+Kσ, k, K = range, and step defines granularity.
	if train == true, use groundtruth dataframe to find ROIs around emissions (e.g. for ER)
"""
function parseFixedFrame(path, ROI_PX, PX_NM, FRAME_NR, groundtruth_dataframe=nothing, range=(1,3), train=false, fixed=0, step=0.1)
    img = Images.load(path)
    NP = size(img)[1]
    @assert(NP > 1)
    lm, lM, μ, σ = peaks(img);
    NegBox, d, score, CC, weightedcentroid, ROIS, modifiedcenter, NegBox, minfiltered, maxfiltered = [nothing for _ in 1:10]
    # gt is a dataframe
    locs_px =nothing
    if groundtruth_dataframe != nothing && train == true
        frame_i = groundtruth_dataframe[groundtruth_dataframe.frame .== FRAME_NR,:]
        locs_px = floor.(Int, [frame_i.xnano frame_i.ynano] ./ PX_NM) .+ 1;
    end
    if fixed == 0
		@assert range[1]<range[2]
		@assert 0 < step < (range[2]-range[1])
		T = range[1]
        for t in range[1]:step:range[2] # For k (paper)
            maxfiltered = [i for i in lM if img[i] >= μ + t*σ]
            minfiltered = [i for i in lm if img[i] <  μ  ]
            d, score, CC, weightedcentroid = findCentroids(ROI_PX, maxfiltered, img, locs_px);
            ROIS, modifiedcenter = findROIS(ROI_PX, weightedcentroid, NP);
            NegBox = findNegatives(ROI_PX, NP, minfiltered, modifiedcenter, locs_px);
            if NegBox != nothing
                break
				T = t
            end
        end
        @assert(NegBox != nothing)
        return img, d, score, CC, weightedcentroid, ROIS, modifiedcenter, NegBox, minfiltered, maxfiltered, T
    else
        maxfiltered = [i for i in lM if img[i] >= fixed]
        return maxfiltered
    end
end

"""
	Computes training data for ER network
	results is a sequence where results[i] matches results (ROIs + metadata) for frame i.
	Writes results to outdir, 1 tiff per ROI, 1 tiff per frame, and 1 tiff per negative ROI, and dataframe to csv.
	PX_NM : 1 pixel is how many nm (assumes isotropic, does not support yet the case where 1 pixel in y != x nm)
"""
function writeTrainingData(results, groundtruthfile, outdir, ROI_PX, PX_NM)
    df = CSV.File(groundtruthfile) |> DataFrames.DataFrame;
    countviolations, cas = Dict(), []
    data = Dict()
    data["metadata"] = ["x nm","y nm","z nm","A photons", "ground truth index"]
    data["data"] = Dict()
    gt=data["data"]
    maxfnr = 0
    @showprogress for (frame_nr, result) in results
        frame_i = df[df.frame .== frame_nr,:]
        GT = [frame_i.xnano frame_i.ynano frame_i.znano frame_i[!, 6] frame_i[!, 1] ]
        frames = [[] for _ in 1:size(GT,1)]
        img, d, score, CC, weightedcentroid, ROIS, modifiedcenter, NegBox, minfiltered, maxfiltered, t = result
        gt[frame_nr] = Dict()
        # Save the ROIs
        for (roindex, roi) in enumerate(ROIS)
            gt[frame_nr][roindex] = []
            xs, ys = roi[:,1], roi[:,2]
            minx, maxx, miny, maxy = minimum(xs),maximum(xs),minimum(ys),maximum(ys)
            @assert maxx - minx == ROI_PX*2 == maxy - miny
            roi_img = img[miny:miny+ROI_PX*2, minx:minx+2*ROI_PX]
            Images.save(joinpath(outdir, "frame_$(frame_nr)_$roindex.tiff"), roi_img)
            gt[frame_nr][roindex] = Dict("roi"=>[], "roi_coord"=>[minx, miny, ROI_PX])
            for row in 1:size(GT,1)
                x, y, z, A, gti = GT[row,:]
                xpx, ypx = floor.(Int, [x y] ./ PX_NM) .+ 1;
                if all((minx, miny) .<= (xpx, ypx) .<= (maxx, maxy))
                    push!(gt[frame_nr][roindex]["roi"], [x, y, z, A, gti])
                    push!(frames[row], roindex)
                end
            end
        end
        # Save the negative ROI
        nroi, roindex=NegBox[1], -1
        xs, ys = [nroi[1][1] nroi[2][1]] , [nroi[1][2] nroi[2][2]]
        minx, maxx, miny, maxy = minimum(xs), maximum(xs), minimum(ys), maximum(ys)
        @assert maxx - minx == maxy - miny
        @assert maxx - minx == ROI_PX*2
        nroi_img = img[miny:miny+ROI_PX*2, minx:minx+2*ROI_PX]
        gt[frame_nr][roindex] = Dict("roi"=>[], "roi_coord"=>[minx, miny, ROI_PX])
        for row in 1:size(GT,1)
            x, y, z, A, gti = GT[row,:]
            xpx, ypx = floor.(Int, [x y] ./ PX_NM).+1;
            if all((minx, miny) .<= (xpx, ypx) .<= (maxx, maxy))
                push!(gt[frame_nr][roindex]["roi"], [x, y, z, A, gti])
                countviolations[frame_nr] = true
                push!(cas, A)
                push!(frames[row], roindex)
            end
        end
        gt[frame_nr]["missed"] = frames
        Images.save(joinpath(outdir, "frame_$(frame_nr).tiff"), img)
        Images.save(joinpath(outdir, "frame_$(frame_nr)_negative.tiff"), nroi_img)
        maxfnr = max(maxfnr, frame_nr)
    end
    data["framecount"] = maxfnr
    open(joinpath(outdir, "trainingdata.json"),"w") do f
        write(f, JSON.json(data));
    end
    return data
end

"""
    Parse a frame into ROIs ROI_PX * 2 + 1 squared
"""
function parseFrame(path, ROI_PX, PX_NM, FRAME_NR, groundtruth_dataframe=nothing, range=(1,3), train=false)
    img = Images.load(path)
    sx, sy = size(img)
    if sx != sy
        @warn "WARNING, non square image $sx x $sy"
        img = img[1:min(sx,sy), 1:min(sx,sy)]
        @warn "TRUNCATED TO $(size(img))"
    end
    NP = size(img,1)
    @assert(NP > 1)
    lm, lM, μ, σ = peaks(img);
    NegBox, d, score, CC, weightedcentroid, ROIS, modifiedcenter, NegBox, minfiltered, maxfiltered = [nothing for _ in 1:10]
    # gt is a dataframe
    locs_px =nothing
    if groundtruth_dataframe != nothing && train == true
        frame_i = groundtruth_dataframe[groundtruth_dataframe.frame .== FRAME_NR,:]
        locs_px = floor.(Int, [frame_i.xnano frame_i.ynano] ./ PX_NM) .+ 1;
    end
    T = 0.0
    for t in range[1]:0.1:range[2]
        T = μ + t*σ
        maxfiltered = [i for i in lM if img[i] >= T]
        minfiltered = [i for i in lm if img[i] <  μ  ]
        d, score, CC, weightedcentroid = findCentroids(ROI_PX, maxfiltered, img, locs_px);
        ROIS, modifiedcenter = findROIS(ROI_PX, weightedcentroid, NP);
        NegBox = findNegatives(ROI_PX, NP, minfiltered, modifiedcenter, locs_px);
        # T = μ + t*σ
        if NegBox != nothing
            break
        end
    end
    @assert(NegBox != nothing)
    return img, d, score, CC, weightedcentroid, ROIS, modifiedcenter, NegBox, minfiltered, maxfiltered, T
end


"""
    findCentrois(radius, maxima, img[, gt=nothing])

    For a given image, radius and maxima finds an optimal covering with tiles radius² px.

    The optimal distance used, the corresponding score, connected components and intensity weighted centroid is returned.

"""
function findCentroids(radius, maxima, img, gt=nothing)
    ROI_PX = radius
    distances = [(ROI_PX/2),(ROI_PX/2)*√(2), ROI_PX, ROI_PX*√(2)]
    vertices = [[i[2], i[1]] for i in maxima];
    if gt != nothing
        vertices = vcat(vertices, [[gt[i,1],gt[i,2]] for i in 1:size(gt,1)])
    end
    adjmat = buildAdjMat(vertices)
    @assert(size(adjmat)[1] == size(vertices)[1])
    results = Dict()
    smax = -Inf
    d = 0.0 # was nothing
    for distance in distances
        graph = buildGraph(adjmat, distance)
        CC = connectedComponents(size(vertices)[1], graph)
        weightedcentroid = computeCentroid(CC, img, vertices)
        score = scoreCC(weightedcentroid, CC, ROI_PX, vertices)
        results[distance] = [CC, weightedcentroid, score]
        smax, d = (score > smax) ? (score , distance) : (smax, d)
    end
    @assert(smax != -Inf)
    CC, weightedcentroid, score = results[d]
    return d, score, CC, weightedcentroid
end

"""
    From a set of points (for which norm(u, v) >= 0), construct the lower triangular matrix representing distance (symmetric) between u,v
"""
function buildAdjMat(vertices)
    N = size(vertices)[1]
    adjmat = Array{Float64}(undef, N, N)
    for i in 1:N, j in i+1:N
        adjmat[i,j] = norm(vertices[j] - vertices[i])
    end
    return adjmat
end

"""
    Scores a tiling for a connected component
    Score >= 0 if tiling is optimal (complete), else < 0
    Returns the mean score of all cc's
"""
function scoreCC(centroids, CC, distance, vertices)
    N = size(CC,1)
    score = 0.0
    A1 = 2*distance*2 # ROI area
    for (centroid, cc) in zip(centroids, CC)
        if size(cc)[1] < 2
            continue
        end
        points_in_cc = vertices[cc]
        mx, my = typemax(Int), typemax(Int)
        Mx, My = typemin(Int), typemin(Int)
        for p in points_in_cc
            x,y = p
            mx = min(x, mx)
            my = min(y, my)
            Mx = max(x, Mx)
            My = max(y, My)
        end
        A2 = max((Mx-mx),1) * max((My-my),1) # bounding box area
        score += A1-A2 # ROI - bounding box, negative if bounding box is larger
        @assert(A2 > 0)
    end
    return score/N
end

"""
	Center the values of array xs s.t. y[i] = (x[i] - minx)/(maxx - minx)
	Note that if max ~ min, even if they are not the same, dividing by a very small floating point number is very unstable in terms of precision.
"""
function centernorm(xs)
	m, M = minimum(xs[:]), maximum(xs[:])
	@assert M != m
	ys = copy(xs)
	ys = (xs .- m) / (M-m)
	return ys
end

"""
	Converts a set of centroids (cs) to ROIs of area (2*ROI_PX +1)^2, taking care to shift the ROI if the center is too close to the edge (no padding is applied)
"""
function ctoroi(cs, ROI_PX)
    ROIS = []
    for wc in cs
		# Convert this line to something less evil
        x, y = xc, yc = wc[2], wc[1]
        xmin, ymin, xmax, ymax = x - ROI_PX, y - ROI_PX, x + ROI_PX, y + ROI_PX
        @assert(ymax - ymin == ROI_PX*2)
        @assert(xmax - xmin == ROI_PX*2)
        xmin, xmax, xc = (xmin < 1) ? (1, x+ROI_PX  - (xmin - 1), x - (xmin - 1)) : (xmin, xmax, x)
        ymin, ymax, yc = (ymin < 1) ? (1, y+ROI_PX  - (ymin - 1), y - (ymin -1)) : (ymin, ymax, y)
        @assert(ymax - ymin == 2 * ROI_PX)
        @assert(xmax - xmin == 2 * ROI_PX)
        ROI = [xmin ymin; xmax ymin; xmin ymax; xmax ymax]
        push!(ROIS, ROI)
    end
    return ROIS
end

"""
    Given adjacency matrix of undirected graph, compute a dictionary based graph with edges <= maxdistance
	Returns a dictionary representation (e.g. sparse), if edge[i,j] > maxdistance then i,j not in graph.
"""
function buildGraph(adjmat, maxdistance)
    graph = Dict()
    N = size(adjmat)[1]
    for i in 1:N, j in i+1:N
        dij = adjmat[i, j]
        if dij <= maxdistance
            get!(graph, i, Dict())
            graph[i][j] = dij
            get!(graph, j, Dict())
            graph[j][i] = dij
        end
    end
    return graph
end


"""
	DFS based CC algorithm
	For |V| large this can trigger stack overflow due to recursion (I hoped LLVM would tail recurse, but not yet)
"""
function connectedComponents(N, adjmat)
    CC = []
    visited = [false for i in 1:N];
    for i in 1:N
        if !visited[i]
            cc = []
            dfs(adjmat, i, visited, cc)
            push!(CC, cc)
        end
    end
    return CC
end

function ifc(images, randomize=false, seed=0, step=1)
    """
    Calculate Pearson inter frame correlation[^1].
    Returns Array{Float64}(N-1) of IFC, means and stds of intensity.
    ...
    # Arguments
    - `images::Collection{String}`: A sequence of filenames to images
    - `randomize::Boolean=false`: Shuffle the sequence (often used to test ifc w.r.t. sequential input)
    - `seed::Integer=0`:seed of the RNG if randomize = true
    - `step::Integer=1`:step size if randomize = false, 1 means ifc(first, second), k is ifc(imgi, imgi+k)
    ...

    [^1]: https://hal.archives-ouvertes.fr/file/index/docid/860912/filename/Image_Processing_Using_Pearsona_s_Correlation_Coefficient_-_Applications_on_Autonomous_Robotics_final.pdf
    """
    N = size(images, 1)
    @assert(N > 0)
    indices = 1:N
    if randomize
        Random.seed!(seed)
        indices = Random.shuffle(indices)
    end
    @assert(N > 0)
    means = Array{Float64}(undef, N)
    stds = Array{Float64}(undef, N)
    imgi, imgj = nothing, nothing
    rijs = Array{Float64}(undef, N-1)
    dist = 0
    for _ith in 1:N-1
        ith = indices[_ith]
        jth = indices[_ith+1]
        dist += abs(ith-jth)
        if imgi == nothing
            imgi = Images.load(images[ith])
            means[_ith] = Statistics.mean(imgi)
            stds[_ith] = Statistics.std(imgi)
        end
        imgj = Images.load(images[jth])
        mean_i = means[_ith]
        means[_ith+1] = mean_j = Statistics.mean(imgj)
        stds[_ith+1] = Statistics.std(imgj)
        rij = sum((imgi .- mean_i).*(imgj .- mean_j)) / (√( sum((imgi .- mean_i).^2)) * √( sum((imgj .- mean_j).^2)) )
        @assert(-1 <= rij <= 1)
        rijs[_ith] = rij
        imgi = imgj
    end
    dist /= (N-1)
    return rijs, means, stds, dist
end

"""
	Shorthand boundary check
"""
function check(coord, m, M, b)
    return m+b < coord < M-b
end

"""
	Find negative (noise) ROIs in a frame
	If gt != nothing, then use ground truth to check our candidates.
	Returns the first ROI we think is noise, if none are found returns nothing (can happen in very high density frames)
	minfiltered is an array of local minima (exceeding a threshold)
	centers are the centroids of already found ROIs
	It's not that difficult to find 'noise' in a frame, but we're trying to find ROIs (e.g. based on local maxima) that are more likely to be noise,
	by checking if they contain a local minima. An ideal PSF is smooth, so would not have a local minima within range of a local maxima (not always true)
	Naturally in high SNR conditions this condition is less and less reliable, but is still worth pursuing.
	If we have ground truth we ensure that any candidate we find really does not contain a source (no matter how weak).
	It's still possible for this function to fail to find any 'negative' ROI even given ground truth, e.g. if the heuristic of local min not near local max fails.
"""
function findNegatives(ROI_PX, NP, minfiltered, centers, gt)
    for (i,mf) in enumerate(minfiltered)
        check = [norm([mf[2] mf[1]]-[j[1] j[2]]) >= √(2) * ROI_PX*2 for j in centers]
        bord = [ROI_PX < mf[2] < NP-ROI_PX ROI_PX < mf[1] < NP-ROI_PX]
        if reduce(&, check, init=true) && reduce(&, bord, init=true) # If we find a negative box ROI_PX * 2 +1 centered on a local minima, that does not overlap with the existing ROIS, and is inside our frame
            candidate = true
            # checked_gt = [abs(mf[2] - gt[gt_index,:][1]) <= ROI_PX && abs(mf[1] - gt[gt_index,:][2]) <= ROI_PX for gt_index in 1:size(gt,1)]
            if gt != nothing # If we have gt, check our negative frame does not contain gt
                for gt_index in 1:size(gt,1)
                    x, y = gt[gt_index,:]
                    if abs(mf[2] - x) <= ROI_PX && abs(mf[1] - y) <= ROI_PX
                        candidate = false
                        break
                    end
                end
            end
            if candidate
                return [[[mf[2]-ROI_PX mf[1]+ROI_PX],[mf[2]+ROI_PX mf[1]-ROI_PX]]]
            end
        end
    end
    return nothing
end

"""
	Recursive depth first search. Can overrun stack, no tail recursion (yet).
"""
function dfs(edges, i, visited, cc)
    visited[i] = true
    push!(cc, i)
    if haskey(edges, i)
        for (target, distance) in edges[i]
            if !visited[target]
                dfs(edges, target, visited, cc)
            end
        end
    end
end

function mcc(TP, TN, FP, FN)
    return (TP*TN - FP*FN)/√((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
end

"""
	For a set of connected components (of centers) in an image, adjust the centroid by weighting the intensity.
	In SMLM it's possible to have multiple weak local maxima (from emissions k frames past, or with signal degraded due to bleeding), and 1 bright local maxima.
	In that case (and if the ROI is not large enough), we want to capture our 'best bet', so shift the centroid towards the intensity weighted center.
	Note: this function can be optimized a lot more than it is now, with vectorized operations, eachindex, CartesianIndex etc, not to mention parallelization.
	Prellocation can also help.
"""
function computeCentroid(CC, img, centers)
    # Centers are xy points of local intensity maxima
    weightedcentroid = []
    for cc in CC
        # For each connected component find an intensity weighted centroid
        wsx, wsy, intens_total = zeros(Float64, 3)
        for point_index in cc
            x, y = centers[point_index]
            intens_total += Float64(img[y,x])
        end
        for point_index in cc
            x, y = centers[point_index]
            intens = Float64(img[y,x]) / intens_total
            wsx += Float64(x) * intens
            wsy += Float64(y) * intens
        end
        push!(weightedcentroid,[round(Int,wsx), round(Int,wsy)])
    end
    return weightedcentroid
end

"""
	Given a set of intensity weighted centroids, fit ROIs around them.
	Simple, but deals with a lot of edge cases/borders etc.
	NP is
"""
function findROIS(ROI_PX, weightedcentroids, NP)
    ROIS = []
    modifiedcenters = []
    for wc in weightedcentroids
        x, y = xc, yc = wc[1], wc[2]
        xmin, ymin, xmax, ymax = x - ROI_PX, y - ROI_PX, x + ROI_PX, y + ROI_PX
        @assert(ymax - ymin == ROI_PX*2)
        @assert(xmax - xmin == ROI_PX*2)
        xmin, xmax, xc = (xmin < 1) ? (1, x+ROI_PX  - (xmin - 1), x - (xmin - 1)) : (xmin, xmax, x)
        xmin, xmax, xc = (xmax > NP) ? (xmin - (xmax - NP), NP, x - (xmax-NP)) : (xmin, xmax, x)
        ymin, ymax, yc = (ymin < 1) ? (1, y+ROI_PX  - (ymin - 1), y - (ymin -1)) : (ymin, ymax, y)
        ymin, ymax, yc = (ymax > NP) ? (ymin - (ymax - NP), NP, y- (ymax-NP)) : (ymin, ymax, y)
        @assert(ymax - ymin == 2 * ROI_PX)
        @assert(xmax - xmin == 2 * ROI_PX)
        push!(modifiedcenters, [xc,yc])
        ROI = [xmin ymin; xmax ymin; xmin ymax; xmax ymax]
        push!(ROIS, ROI)
    end
    return ROIS, modifiedcenters
end

#That's all folks
end # module
