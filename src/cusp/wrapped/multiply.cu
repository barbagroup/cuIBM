#include <cusp/wrapped/multiply.h>

namespace cusp
{
namespace wrapped
{
  // SpMM
  template <class Matrix, class Vector>
  void multiply(const Matrix &A, const Vector &B, Vector &C)
  {
    cusp::multiply(A,B,C);
  }

// Explicit instantiation
// Matrix-Matrix
template void multiply<cooD,cooD>(const cooD &A, const cooD &B, cooD &C);
template void multiply<csrD,csrD>(const csrD &A, const csrD &B, csrD &C);
template void multiply<cooH,cooH>(const cooH &A, const cooH &B, cooH &C);
template void multiply<csrH,csrH>(const csrH &A, const csrH &B, csrH &C);

// Matrix-Vector
template void multiply<cooD,vecD>(const cooD &A, const vecD &B, vecD &C);
template void multiply<csrD,vecD>(const csrD &A, const vecD &B, vecD &C);
template void multiply<cooH,vecH>(const cooH &A, const vecH &B, vecH &C);
template void multiply<csrH,vecH>(const csrH &A, const vecH &B, vecH &C);

} // end namespace wrapped
} // end namespace cusp
