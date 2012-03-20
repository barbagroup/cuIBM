#include <cusp/wrapped/subtract.h>

namespace cusp
{
namespace wrapped
{
  template <class Matrix>
  void subtract(const Matrix &A, const Matrix &B, Matrix &C)
  {
    cusp::subtract(A,B,C);
  }

  template void subtract<cooD>(const cooD &A, const cooD &B, cooD &C);
  template void subtract<csrD>(const csrD &A, const csrD &B, csrD &C);
  template void subtract<cooH>(const cooH &A, const cooH &B, cooH &C);
  template void subtract<csrH>(const csrH &A, const csrH &B, csrH &C);
} // end namespace wrapped
} // end namespace cusp
