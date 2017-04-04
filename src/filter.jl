include("variance.jl")
import ImageFilters

function getamean(variancenorm::Array{Float64,3}, covariance::Array{Float64,3},
    boxradius::Int64)
  channels, imagewidth, imageheight = size(covariance)
  amean = Array{Float64}(3, imagewidth, imageheight)
  for pixely in 1:imageheight
    for pixelx  in 1:imagewidth
      @inbounds redcov = covariance[1, pixelx, pixely]
      @inbounds greencov = covariance[2, pixelx, pixely]
      @inbounds bluecov = covariance[3, pixelx, pixely]
      amean[1, pixelx, pixely] = (redcov * variancenorm[1, pixelx, pixely] +
        greencov * variancenorm[2, pixelx, pixely] + bluecov * variancenorm[3,
        pixelx, pixely])
      amean[2, pixelx, pixely] = (redcov * variancenorm[2, pixelx, pixely] +
        greencov * variancenorm[4, pixelx, pixely] + bluecov * variancenorm[5,
        pixelx, pixely])
      amean[3, pixelx, pixely] = (redcov * variancenorm[3, pixelx, pixely] +
        greencov * variancenorm[5, pixelx, pixely] + bluecov * variancenorm[6,
        pixelx, pixely])
    end
  end
  return amean
end

function computemeandiff!(amean::Array{Float64,3}, bmean::Array{Float64,2},
    guidancemean::Array{Float64,3})
  channels, imagewidth, imageheight = size(guidancemean)
  for pixely in 1:imageheight
    for pixelx in 1:imagewidth
      @inbounds bmean[pixelx, pixely] -= (amean[1, pixelx, pixely] *
        guidancemean[1, pixelx, pixely] + amean[2, pixelx, pixely] *
        guidancemean[2, pixelx, pixely] + amean[3, pixelx, pixely] *
        guidancemean[3, pixelx, pixely])
    end
  end
end

function getfilteredimage(amean::Array{Float64,3}, bmean::Array{Float64,2},
    guidanceimage::Array{Float64,3}, boxradius::Int64)
  channels, imagewidth, imageheight = size(guidanceimage)
  filteredimage = Array{Float64}(imagewidth, imageheight)

  for pixely in 1:imageheight
    for pixelx in 1:imagewidth
      @inbounds filteredimage[pixelx, pixely] = ((amean[1, pixelx, pixely] *
        guidanceimage[1, pixelx, pixely]) + (amean[2, pixelx, pixely] *
        guidanceimage[2, pixelx, pixely]) + (amean[3, pixelx, pixely] *
        guidanceimage[3, pixelx, pixely]) + bmean[pixelx, pixely])
    end
  end
  return filteredimage
end

function getguidedfilter(guidanceimage::Array{Float64,3},
    guidancemean::Array{Float64,3}, variancenorm::Array{Float64,3},
    targetimage::Array{Float64,2}, boxradius::Int64, epsilon::Float64)
  targetmean = ImageFilters.boxfilter(targetimage, boxradius)
  covariance = getcovariance(guidanceimage, targetimage, guidancemean,
    targetmean, boxradius)
  amean = getamean(variancenorm, covariance, boxradius)
  bmean = targetmean
  computemeandiff!(amean, bmean, guidancemean)
  amean = ImageFilters.boxfilter(amean, boxradius)
  bmean = ImageFilters.boxfilter(bmean, boxradius)
  return getfilteredimage(amean, bmean, guidanceimage, boxradius)
end

function getguidedfilter(guidanceimage::Array{Float64, 3},
    targetimage::Array{Float64, 2}, boxradius::Int64, epsilon::Float64)
  guidancemean, variancenorm = initiatefilter(guidanceimage, boxradius, epsilon)
  covariance = getcovariance(guidanceimage, targetimage, guidancemean,
    targetmean, boxradius)
  amean = getamean(variancenorm, covariance, boxradius)
  bmean = targetmean
  computemeandiff!(amean, bmean, guidancemean)
  amean = ImageFilters.boxfilter(amean, boxradius)
  bmean = ImageFilters.boxfilter(bmean, boxradius)
  return getfilteredimage(amean, bmean, guidanceimage, boxradius)
end
