function getmatrixinverse!(matrix::Array{Float32,2},
    inversematrix::Array{Float32,2})
  @inbounds a11 = ((matrix[3,3]*matrix[2,2]) - (matrix[2,3]*matrix[2,3]))
  @inbounds a12 = ((matrix[1,3]*matrix[2,3]) - (matrix[3,3]*matrix[1,2]))
  @inbounds a13 = ((matrix[1,2]*matrix[2,3]) - (matrix[1,3]*matrix[2,2]))
  @inbounds determinant = ((matrix[1,1]*a11) + (matrix[2,1]*a12) +
    (matrix[1,3]*a13))

  @inbounds inversematrix[1,1] = a11 / determinant
  @inbounds inversematrix[1,2] = a12 / determinant
  @inbounds inversematrix[1,3] = a13 / determinant
  @inbounds inversematrix[2,2] = ((matrix[3,3]*matrix[1,1]) -
    (matrix[1,3]*matrix[1,3])) / determinant
  @inbounds inversematrix[2,3] = ((matrix[1,2]*matrix[1,3]) -
    (matrix[1,1]*matrix[2,3])) / determinant
  @inbounds inversematrix[3,3] = ((matrix[1,1]*matrix[2,2]) -
  (matrix[1,2]*matrix[1,2])) / determinant
  @inbounds inversematrix[2,1] = a12 / determinant
  @inbounds inversematrix[3,1] = a13 / determinant
  @inbounds inversematrix[3,2] = inversematrix[2,3]
end
