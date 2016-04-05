/***************************************************************************//**
 * \file preconditioner.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c preconditioner.
 */


#pragma once

#include <cusp/linear_operator.h>
#include <cusp/csr_matrix.h>
#include <cusp/precond/diagonal.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/precond/ainv.h>
#include <cusp/format_utils.h>
#include "types.h"


/**
 * \class preconditioner
 * \brief Stores the preconditioner for a given system.
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

	// constructor
	preconditioner();

	// constructor overloading
	preconditioner(const Matrix &A, preconditionerType _type);

	// destructor
	~preconditioner();

	// update the preconditioner of the system
	void update(const Matrix &A);

	// overload the operator ()
	template <typename VectorType1, typename VectorType2>
	void operator()(const VectorType1 &x, VectorType2 &y) const;
};

/**
 * \brief Constructor.
 */
template <class Matrix>
preconditioner<Matrix>::preconditioner()
{
}

/**
 * \brief Constructor overloading. Computes the preconditioner 
 *        given a matrix and a type of preconditioner.
 *
 * \param A matrix of the system (instance of the class \c Matrix)
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
	if (type == AINV)
		LO = new cusp::precond::nonsym_bridson_ainv<ValueType, MemorySpace>(A);
	else
	{
		std::cout << "ERROR: Choose a valid preconditioner!\n" << std::endl;
		exit(0);
	}
}

/**
 * \brief Destructor. Deletes the preocnditioner.
 */
template <typename Matrix>
preconditioner<Matrix>::~preconditioner()
{
	delete LO;
}

/**
 * \brief Updates the preconditioner of the system.
 *
 * \param A matrix of the system (instance of the class \c Matrix)
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
	if (type == AINV)
		*LO = cusp::precond::nonsym_bridson_ainv<ValueType, MemorySpace>(A);
	else
	{
		std::cout << "ERROR: Choose a valid preconditioner!\n" << std::endl;
		exit(0);
	}
}

/**
 * \brief Overloads the operator (). This is required due to the way preconditioners are 
          implemented in Cusp - as linear operators on a vector.
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
	else if (type == AINV)
	{
		cusp::precond::nonsym_bridson_ainv<value_type, memory_space> *AI =
			static_cast<cusp::precond::nonsym_bridson_ainv<value_type, memory_space> *>(LO);
		AI->operator()(x,y);
	}
	else
	{
		printf("ERROR: Choose a valid preconditioner!\n");
		exit(0);
	}
}
