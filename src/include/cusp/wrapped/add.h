#pragma once

#include <types.h>
#include <cusp/elementwise.h>

namespace cusp
{
namespace wrapped
{
  template <class Matrix>
  void add(const Matrix &A, const Matrix &B, Matrix &C);

} // end namespace wrapped
} // end namespace cusp
