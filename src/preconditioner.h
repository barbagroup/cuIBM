/***************************************************************************//**
* \file  preconditioner.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration and definition of the class \c preconditioner
*/

#pragma once

#include "types.h"
#include <cusp/linear_operator.h>
#include <cusp/csr_matrix.h>
#include <cusp/precond/diagonal.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/format.h>

/**
* \class preconditioner
* \brief Store the preconditioner for a given system
*/
template <typename Matrix>
class preconditioner
{
	preconditionerType type;	///< type of preconditioner

	cusp::linear_operator<typename Matrix::value_type,
	                      typename Matrix::memory_space,
	                      typename Matrix::index_type>* LO;

public:
	typedef typename Matrix::index_type   index_type;
	typedef typename Matrix::value_type   value_type;
	typedef typename Matrix::memory_space memory_space;
	typedef typename cusp::unknown_format format;

	/**
	* \brief Constructor of the class \c preconditioner
	*/
	preconditioner();

	/**
	* \brief Overloaded constructor of the class \c preconditioner
	*/
	preconditioner(const Matrix &A, preconditionerType _type);

	/**
	* \brief Destructor of the class \c preconditioner
	*/
	~preconditioner();

	/**
	* \brief Update the preconditioner of the system
	*/
	void update(const Matrix &A);

	/**
	* \brief Overload the operator ()
	*/
	template <typename VectorType1, typename VectorType2>
	void operator()(const VectorType1 &x, VectorType2 &y) const;
};

/**
* \brief Constructor of the class \c preconditioner
*/
template <class Matrix>
preconditioner<Matrix>::preconditioner()
{
}

/**
* \brief Overloaded constructor of the class \c preconditioner
*
* \param A matrix of the system
* \param _type the type of preconditioner
*/
template <class Matrix>
preconditioner<Matrix>::preconditioner(const Matrix &A, preconditionerType _type)
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
		LO = new cusp::precond::aggregation::smoothed_aggregation<IndexType, ValueType, MemorySpace>(A);
	else
	{
		std::cout << "ERROR: Choose a valid preconditioner!\n" << std::endl;
		exit(0);
	}
}

/**
* \brief Destructor of the class \c preconditioner
*/
template <typename Matrix>
preconditioner<Matrix>::~preconditioner()
{
	delete LO;
}

/**
* \brief Update the preconditioner of the system
*
* \param A matrix of the system
*/
template <class Matrix>
void preconditioner<Matrix>::update(const Matrix &A)
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
		*LO = cusp::precond::aggregation::smoothed_aggregation<IndexType, ValueType, MemorySpace>(A);
	else
	{
		std::cout << "ERROR: Choose a valid preconditioner!\n" << std::endl;
		exit(0);
	}
}

// Why is this required? Need to look into the implementation of preconditioners in Cusp.
/**
* \brief Overload of the operator ()
*/
template <typename Matrix>
template <typename VectorType1, typename VectorType2>
void preconditioner<Matrix>::operator()(const VectorType1 &x, VectorType2 &y) const
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
		cusp::precond::aggregation::smoothed_aggregation<index_type, value_type, memory_space> *SA =
			static_cast<cusp::precond::aggregation::smoothed_aggregation<index_type, value_type, memory_space> *>(LO);
		SA->operator()(x,y);
	}
	else
	{
		printf("ERROR: Choose a valid preconditioner!\n");
		exit(0);
	}
}