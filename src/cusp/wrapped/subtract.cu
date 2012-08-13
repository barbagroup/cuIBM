/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
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
