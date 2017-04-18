import ImageFilters
include("matrixutils.jl")

function getvariance(guidanceimage::Array{Float32,3},
    guidancemean::Array{Float32,3}, boxradius::Int64)
  channels, imagewidth, imageheight = size(guidanceimage)
  greenbluevariance = Array{Float32}(3, imagewidth, imageheight)
  redvariance = Array{Float32}(3, imagewidth, imageheight)
  for pixely in 1:imageheight
    for pixelx in 1:imagewidth
      redvariance[1, pixelx, pixely] = (guidanceimage[1, pixelx, pixely] *
        guidanceimage[1, pixelx, pixely])
      redvariance[2, pixelx, pixely] = (guidanceimage[1, pixelx, pixely] *
        guidanceimage[2, pixelx, pixely])
      redvariance[3, pixelx, pixely] = (guidanceimage[1, pixelx, pixely] *
        guidanceimage[3, pixelx, pixely])
      greenbluevariance[1, pixelx, pixely] = (guidanceimage[2, pixelx, pixely] *
        guidanceimage[2, pixelx, pixely])
      greenbluevariance[2, pixelx, pixely] = (guidanceimage[2, pixelx, pixely] *
        guidanceimage[3, pixelx, pixely])
      greenbluevariance[3, pixelx, pixely] = (guidanceimage[3, pixelx, pixely] *
        guidanceimage[3, pixelx, pixely])
    end
  end
  redvariance = ImageFilters.boxfilter(redvariance, (boxradius, boxradius))
  greenbluevariance = ImageFilters.boxfilter(greenbluevariance, (boxradius,
    boxradius))

  for pixely in 1:imageheight
    for pixelx in 1:imagewidth
      redvariance[1, pixelx, pixely] -= (guidancemean[1, pixelx, pixely] *
        guidancemean[1, pixelx, pixely])
      redvariance[2, pixelx, pixely] -= (guidancemean[1, pixelx, pixely] *
        guidancemean[2, pixelx, pixely])
      redvariance[3, pixelx, pixely] -= (guidancemean[1, pixelx, pixely] *
        guidancemean[3, pixelx, pixely])
      greenbluevariance[1, pixelx, pixely] -= (guidancemean[2, pixelx, pixely]
        * guidancemean[2, pixelx, pixely])
      greenbluevariance[2, pixelx, pixely] -= (guidancemean[2, pixelx, pixely]
        * guidancemean[3, pixelx, pixely])
      greenbluevariance[3, pixelx, pixely] -= (guidancemean[3, pixelx, pixely]
        * guidancemean[3, pixelx, pixely])
    end
  end
  return redvariance, greenbluevariance
end

function getvariancenorm(redvariance::Array{Float32,3},
    greenbluevariance::Array{Float32,3} ,epsilon::Float32)
  channels, imagewidth, imageheight = size(redvariance)
  variancenorm = Array{Float32}(6, imagewidth, imageheight)
  sigma = Array{Float32}(3,3)
  invsigma = Array{Float32}(3,3)
  for pixely in 1:imageheight
    for pixelx in 1:imagewidth
      @inbounds sigma[1,1] = redvariance[1, pixelx, pixely] + epsilon
      @inbounds sigma[1,2] = redvariance[2, pixelx, pixely]
      @inbounds sigma[1,3] = redvariance[3, pixelx, pixely]
      @inbounds sigma[2,1] = redvariance[2, pixelx, pixely]
      @inbounds sigma[2,2] = greenbluevariance[1, pixelx, pixely] + epsilon
      @inbounds sigma[2,3] = greenbluevariance[2, pixelx, pixely]
      @inbounds sigma[3,1] = redvariance[3, pixelx, pixely]
      @inbounds sigma[3,2] = greenbluevariance[2, pixelx, pixely]
      @inbounds sigma[3,3] = greenbluevariance[3, pixelx, pixely] + epsilon

      getmatrixinverse!(sigma, invsigma)

      @inbounds variancenorm[1, pixelx, pixely] = invsigma[1,1]
      @inbounds variancenorm[2, pixelx, pixely] = invsigma[1,2]
      @inbounds variancenorm[3, pixelx, pixely] = invsigma[1,3]
      @inbounds variancenorm[4, pixelx, pixely] = invsigma[2,2]
      @inbounds variancenorm[5, pixelx, pixely] = invsigma[2,3]
      @inbounds variancenorm[6, pixelx, pixely] = invsigma[3,3]
    end
  end
  return variancenorm
end

function initiatefilter(guidanceimage::Array{Float32,3},
    boxradius::Int64, epsilon::Float32)
  guidancemean = ImageFilters.boxfilter(guidanceimage, (boxradius,boxradius))
  redvariance, greenbluevariance = getvariance(guidanceimage, guidancemean,
    boxradius)
  variancenorm = getvariancenorm(redvariance, greenbluevariance, epsilon)
  return guidancemean, variancenorm
end
