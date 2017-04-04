import ImageFilters

function getcovariance(guidanceimage::Array{Float64,3},
    targetimage::Array{Float64,2}, guidancemean::Array{Float64,3},
    targetmean::Array{Float64,2}, boxradius::Int64)
  channels, imagewidth, imageheight = size(guidancemean)
  covariance = Array{Float64}(channels, imagewidth, imageheight)

  productimage = Array{Float64}(channels, imagewidth, imageheight)
  for pixely in 1:imageheight
    for  pixelx in 1:imagewidth
      @inbounds productimage[1, pixelx, pixely] = (guidanceimage[1, pixelx,
        pixely] * targetimage[pixelx, pixely])
      @inbounds productimage[2, pixelx, pixely] = (guidanceimage[2, pixelx,
        pixely] * targetimage[pixelx, pixely])
      @inbounds productimage[3, pixelx, pixely] = (guidanceimage[3, pixelx,
        pixely] * targetimage[pixelx, pixely])
    end
  end
  productimage = ImageFilters.boxfilter(productimage, boxradius)

  for pixely in 1:imageheight
    for pixelx in 1:imagewidth
      @inbounds covariance[1, pixelx, pixely] = (productimage[1, pixelx, pixely]
        - (guidancemean[1, pixelx, pixely] * targetmean[pixelx, pixely]))
      @inbounds covariance[2, pixelx, pixely] = (productimage[2, pixelx, pixely]
        - (guidancemean[2, pixelx, pixely] * targetmean[pixelx, pixely]))
      @inbounds covariance[3, pixelx, pixely] = (productimage[3, pixelx, pixely]
        - (guidancemean[3, pixelx, pixely] * targetmean[pixelx, pixely]))
    end
  end
  return covariance
end
