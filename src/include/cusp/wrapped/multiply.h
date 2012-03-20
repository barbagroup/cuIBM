#pragma once

#include <cusp/csr_matrix.h>
#include <cusp/coo_matrix.h>
#include <cusp/multiply.h>
#include <types.h>

namespace cusp
{
namespace wrapped
{
  template <class Matrix, class Vector>
  void multiply(const Matrix &A, const Vector &B, const Vector &C);

} // end namespace wrapped
} // end namespace cusp

