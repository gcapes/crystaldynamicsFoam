/*---------------------------------------------------------------------------*\
  
  FASTOR UTILITY FUNCTIONS FOR OPENFOAM
  ======================================
  
  Conversion between OpenFOAM tensors and Fastor tensors
  
\*---------------------------------------------------------------------------*/

#ifndef FASTOR_OPENFOAM_UTILS_H
#define FASTOR_OPENFOAM_UTILS_H

#include "tensor.H"
#include "Fastor/Fastor.h"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convert OpenFOAM tensor to Fastor 3x3 tensor
inline Fastor::Tensor<double, 3, 3> toFastor(const tensor& T)
{
    Fastor::Tensor<double, 3, 3> F;
    F(0, 0) = T.xx(); F(0, 1) = T.xy(); F(0, 2) = T.xz();
    F(1, 0) = T.yx(); F(1, 1) = T.yy(); F(1, 2) = T.yz();
    F(2, 0) = T.zx(); F(2, 1) = T.zy(); F(2, 2) = T.zz();
    return F;
}

// Convert Fastor 3x3 tensor to OpenFOAM tensor
inline tensor toOpenFOAM(const Fastor::Tensor<double, 3, 3>& F)
{
    return tensor
    (
        F(0, 0), F(0, 1), F(0, 2),
        F(1, 0), F(1, 1), F(1, 2),
        F(2, 0), F(2, 1), F(2, 2)
    );
}

// Convert OpenFOAM vector to Fastor 3-vector
inline Fastor::Tensor<double, 3> toFastor(const vector& v)
{
    Fastor::Tensor<double, 3> F;
    F(0) = v.x();
    F(1) = v.y();
    F(2) = v.z();
    return F;
}

// Convert Fastor 3-vector to OpenFOAM vector
inline vector toOpenFOAM(const Fastor::Tensor<double, 3>& F)
{
    return vector(F(0), F(1), F(2));
}

// Compute trace using Fastor
inline scalar traceFastor(const tensor& T)
{
    auto F = toFastor(T);
    return Fastor::trace(F);
}

// Compute determinant using Fastor
inline scalar detFastor(const tensor& T)
{
    auto F = toFastor(T);
    return Fastor::determinant(F);
}

// Compute inverse using Fastor
inline tensor invFastor(const tensor& T)
{
    auto F = toFastor(T);
    auto Finv = Fastor::inverse(F);
    return toOpenFOAM(Finv);
}

// Compute transpose using Fastor
inline tensor transFastor(const tensor& T)
{
    auto F = toFastor(T);
    auto Ft = Fastor::transpose(F);
    return toOpenFOAM(Ft);
}

// Double contraction: A:B (sum over both indices)
inline scalar doubleContraction
(
    const tensor& A,
    const tensor& B
)
{
    auto FA = toFastor(A);
    auto FB = toFastor(B);
    
    // A:B = sum_ij A_ij * B_ij
    return Fastor::sum(FA * FB);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // FASTOR_OPENFOAM_UTILS_H

// ************************************************************************* //