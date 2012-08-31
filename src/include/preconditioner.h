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

/**
* @file  preconditioner.h
* @brief Stores the preconditioner for a given system.
*/

#pragma once

#include <preconditioner.h>
#include <cusp/linear_operator.h>
#include <cusp/csr_matrix.h>
#include <cusp/precond/diagonal.h>
#include <cusp/precond/smoothed_aggregation.h>
#include <cusp/format.h>

template <typename Matrix, typename Vector>
class preconditioner
{	
	preconditionerType type;
	
	cusp::linear_operator<typename Matrix::value_type,
	                      typename Matrix::memory_space,
	                      typename Matrix::index_type>* LO;

public:
	typedef typename Matrix::index_type   index_type;
	typedef typename Matrix::value_type   value_type;
	typedef typename Matrix::memory_space memory_space;
	typedef typename cusp::unknown_format format;
	
	// constructors
	preconditioner();
	preconditioner(const Matrix &A, preconditionerType _type);
	
	// destructor
	~preconditioner();
	
	void update(const Matrix &A);

	// () operator
	template <typename VectorType1, typename VectorType2>
	void operator()(const VectorType1 &x, VectorType2 &b) const;
};

// constructors
template <class Matrix, class Vector>
preconditioner<Matrix,Vector>::preconditioner()
{
}

// this is simple enough
template <class Matrix, class Vector>
preconditioner<Matrix,Vector>::preconditioner(const Matrix &A, preconditionerType _type)
{
	typedef typename Matrix::value_type   ValueType;
	typedef typename Matrix::index_type   IndexType;
	typedef typename Matrix::memory_space MemorySpace;

	type = _type;

	// generate an instance of linear_operator with the required derived class
	// depending on the second parameter
	/*if (type == NONE)
		LO = new cusp::identity_operator<ValueType, MemorySpace>(A.num_rows, A.num_cols);
	else*/
	if (type == DIAGONAL)
		LO = new cusp::precond::diagonal<ValueType, MemorySpace>(A);
	else
	if (type == SMOOTHED_AGGREGATION)
		LO = new cusp::precond::smoothed_aggregation<IndexType, ValueType, MemorySpace>(A);
	else
	{
		std::cout << "ERROR: Choose a valid preconditioner!\n" << std::endl;
		exit(0);
	}
}

// destructor
template <typename Matrix, typename Vector>
preconditioner<Matrix,Vector>::~preconditioner()
{
	delete LO;
}

template <class Matrix, class Vector>
void preconditioner<Matrix,Vector>::update(const Matrix &A)
{
	typedef typename Matrix::value_type   ValueType;
	typedef typename Matrix::index_type   IndexType;
	typedef typename Matrix::memory_space MemorySpace;

	/*if (type == NONE)
		LO = cusp::identity_operator<ValueType, MemorySpace>(A.num_rows, A.num_cols);
	else*/
	if (type == DIAGONAL)
		*LO = cusp::precond::diagonal<ValueType, MemorySpace>(A);
	else
	if (type == SMOOTHED_AGGREGATION)
		*LO = cusp::precond::smoothed_aggregation<IndexType, ValueType, MemorySpace>(A);
	else
	{
		std::cout << "ERROR: Choose a valid preconditioner!\n" << std::endl;
		exit(0);
	}
}

// the operator defined here is ()
// Why is this required? Need to look into the implementation of preconditioners in Cusp.
template <typename Matrix, typename Vector>
template <typename VectorType1, typename VectorType2> 
void preconditioner<Matrix,Vector>::operator()(const VectorType1 &x, VectorType2 &y) const
{

	/*if (type == NONE)
	{
		cusp::identity_operator<value_type, memory_space> *identity = 
			static_cast<cusp::identity_operator<value_type, memory_space> *>(LO);
		printf("dispatching identity->operator()\n");
		identity->operator()(x,b);
	}
	else*/
	if (type == DIAGONAL)
	{
		cusp::precond::diagonal<value_type, memory_space> *diag = 
			static_cast<cusp::precond::diagonal<value_type, memory_space> *>(LO);
		diag->operator()(x,y);
	}
	else if (type == SMOOTHED_AGGREGATION)
	{
		cusp::precond::smoothed_aggregation<index_type, value_type, memory_space> *SA = 
			static_cast<cusp::precond::smoothed_aggregation<index_type, value_type, memory_space> *>(LO);
		SA->operator()(x,y);
	}
	else
	{
		printf("ERROR: Choose a valid preconditioner!\n");
		exit(0);
	}
}
