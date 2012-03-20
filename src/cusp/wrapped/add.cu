#include <cusp/wrapped/add.h>

namespace cusp
{
namespace wrapped
{
  template <class Matrix>
  void add(const Matrix &A, const Matrix &B, Matrix &C)
  {
    cusp::add(A,B,C);
  }

  template void add<cooD>(const cooD &A, const cooD &B, cooD &C);
  template void add<csrD>(const csrD &A, const csrD &B, csrD &C);
  template void add<cooH>(const cooH &A, const cooH &B, cooH &C);
  template void add<csrH>(const csrH &A, const csrH &B, csrH &C);
} // end namespace wrapped
} // end namespace cusp
