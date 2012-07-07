/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

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
