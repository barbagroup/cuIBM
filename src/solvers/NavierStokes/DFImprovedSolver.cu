#include "DFImprovedSolver.h"
#include <solvers/NavierStokes/kernels/generateQT.h>
#include <cusp/io/matrix_market.h>
#include <cusp/blas/blas.h>

template <typename memoryType>
DFImprovedSolver<memoryType>::DFImprovedSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template <typename memoryType>
void DFImprovedSolver<memoryType>::generateQT()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
	    ny = NavierStokesSolver<memoryType>::domInfo->ny;

	const int N_u = (nx-1)*ny;

	NavierStokesSolver<memoryType>::generateQT();
	DirectForcingSolver<memoryType>::updateQ();

	cusp::coo_matrix<int, real, host_memory> QTHost(nx*ny, (nx-1)*ny+nx*(ny-1), 4*nx*ny-2*(nx+ny));
	cusp::blas::fill(QTHost.row_indices, -1);
	cusp::blas::fill(QTHost.column_indices, -1);
	cusp::blas::fill(QTHost.values, 0.0);

	int idx = 0;
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			int row = j*nx+i;
			if(i>0)
			{
				int I = j*(nx-1)+(i-1);
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					QTHost.row_indices[idx] = row;
					QTHost.column_indices[idx] = I;
					QTHost.values[idx] = 1.0;
					idx++;
				}
				else
				{
					bool flag = false;
					int start;
					start = (idx>4)? idx-4 : 0;
					for(int l=start; l<idx && !flag; l++)
					{
						if(QTHost.row_indices[l]==row && QTHost.column_indices[l]==DirectForcingSolver<memoryType>::tags[I])
						{
							flag = true;
							QTHost.values[l] += DirectForcingSolver<memoryType>::coeffs[I];
						}
					}
					if(!flag)
					{
						QTHost.row_indices[idx]    = row;
						QTHost.column_indices[idx] = DirectForcingSolver<memoryType>::tags[I];
						QTHost.values[idx]         = DirectForcingSolver<memoryType>::coeffs[I];
						idx++;
					}
				}
			}
			if(i<nx-1)
			{
				int I = j*(nx-1)+i;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					QTHost.row_indices[idx] = row;
					QTHost.column_indices[idx] = I;
					QTHost.values[idx] = -1.0;
					idx++;
				}
				else
				{
					bool flag = false;
					int start;
					start = (idx>4)? idx-4 : 0;
					for(int l=start; l<idx && !flag; l++)
					{
						if(QTHost.row_indices[l]==row && QTHost.column_indices[l]==DirectForcingSolver<memoryType>::tags[I])
						{
							flag = true;
							QTHost.values[l] -= DirectForcingSolver<memoryType>::coeffs[I];
						}
					}
					if(!flag)
					{
						QTHost.row_indices[idx]    = row;
						QTHost.column_indices[idx] = DirectForcingSolver<memoryType>::tags[I];
						QTHost.values[idx]         = -DirectForcingSolver<memoryType>::coeffs[I];
						idx++;
					}
				}
			}
			if(j>0)
			{
				int I = (j-1)*nx+i+N_u;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					QTHost.row_indices[idx] = row;
					QTHost.column_indices[idx] = I;
					QTHost.values[idx] = 1.0;
					idx++;
				}
				else
				{
					bool flag = false;
					int start;
					start = (idx>4)? idx-4 : 0;
					for(int l=start; l<idx && !flag; l++)
					{
						if(QTHost.row_indices[l]==row && QTHost.column_indices[l]==DirectForcingSolver<memoryType>::tags[I])
						{
							flag = true;
							QTHost.values[l] += DirectForcingSolver<memoryType>::coeffs[I];
						}
					}
					if(!flag)
					{
						QTHost.row_indices[idx]    = row;
						QTHost.column_indices[idx] = DirectForcingSolver<memoryType>::tags[I];
						QTHost.values[idx]         = DirectForcingSolver<memoryType>::coeffs[I];
						idx++;
					}
				}
			}
			if(j<ny-1)
			{
				int I = j*nx+i+N_u;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					QTHost.row_indices[idx] = row;
					QTHost.column_indices[idx] = I;
					QTHost.values[idx] = -1.0;
					idx++;
				}
				else
				{
					bool flag = false;
					int start;
					start = (idx>4)? idx-4 : 0;
					for(int l=start; l<idx && !flag; l++)
					{
						if(QTHost.row_indices[l]==row && QTHost.column_indices[l]==DirectForcingSolver<memoryType>::tags[I])
						{
							flag = true;
							QTHost.values[l] -= DirectForcingSolver<memoryType>::coeffs[I];
						}
					}
					if(!flag)
					{
						QTHost.row_indices[idx]    = row;
						QTHost.column_indices[idx] = DirectForcingSolver<memoryType>::tags[I];
						QTHost.values[idx]         = -DirectForcingSolver<memoryType>::coeffs[I];
						idx++;
					}
				}
			}
		}
	}
	NavierStokesSolver<memoryType>::QT = QTHost;
	NavierStokesSolver<memoryType>::QT.resize(nx*ny, (nx-1)*ny+nx*(ny-1), idx);
	/*std::cout << "\nQT stuff:\n";
	std::cout << "Copied and resized matrix." << std::endl;
	std::cout << "Original size: " << QTHost.values.size() << std::endl;
	std::cout << "Actual size  : " << idx << std::endl;
	cusp::io::write_matrix_market_file(NavierStokesSolver<memoryType>::Q, "Q.mtx");
	std::cout << "Wrote Q to file." << std::endl;
	cusp::io::write_matrix_market_file(NavierStokesSolver<memoryType>::QT, "QT.mtx");
	std::cout << "Wrote QT to file." << std::endl;*/
}

/*
template <typename memoryType>
void DFImprovedSolver<memoryType>::generateC()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
	    ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;
	real         dt  = db["simulation"]["dt"].get<real>();

	const int ii = 2, jj = 2;
	const int N_u = (nx-1)*ny;
	
	int  isColumnNonZero[5][5];
	for(int m=-2; m<=2; m++)
	{
		for(int l=-2; l<=2; l++)
		{
			isColumnNonZero[jj+m][ii+l]=0;
		}
	}

	int num_nonzeros = 0;
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			if(j>0)
			{
				int I = (j-1)*nx+i+N_u;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj-1][ii] += 1;
					isColumnNonZero[jj][ii]   += 1;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					if(diff < -1)
					{
						isColumnNonZero[jj-2][ii] += 1;
						isColumnNonZero[jj-1][ii] += 1;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj-1][ii-1] += 1;
						isColumnNonZero[jj][ii-1]   += 1;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj-1][ii+1] += 1;
						isColumnNonZero[jj][ii+1]   += 1;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj][ii]   += 1;
						isColumnNonZero[jj+1][ii] += 1;
					}
				}
			}
			if(i>0)
			{
				int I = j*(nx-1)+(i-1);
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj][ii-1] += 1;
					isColumnNonZero[jj][ii]   += 1;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					if(diff < -1)
					{
						isColumnNonZero[jj-1][ii-1] += 1;
						isColumnNonZero[jj-1][ii]   += 1;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj][ii-2] += 1;
						isColumnNonZero[jj][ii-1] += 1;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj][ii]   += 1;
						isColumnNonZero[jj][ii+1] += 1;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj+1][ii-1] += 1;
						isColumnNonZero[jj+1][ii]   += 1;
					}
				}
			}
			if(i<nx-1)
			{
				int I = j*(nx-1)+i;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj][ii+1] += 1;
					isColumnNonZero[jj][ii]   += 1;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					if(diff < -1)
					{
						isColumnNonZero[jj-1][ii+1] += 1;
						isColumnNonZero[jj-1][ii]   += 1;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj][ii] += 1;
						isColumnNonZero[jj][ii-1] += 1;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj][ii+2] += 1;
						isColumnNonZero[jj][ii+1] += 1;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj+1][ii+1] += 1;
						isColumnNonZero[jj+1][ii]   += 1;
					}
				}
			}
			if(j<ny-1)
			{
				int I = j*nx+i+N_u;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj+1][ii] += 1;
					isColumnNonZero[jj][ii]   += 1;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					if(diff < -1)
					{
						isColumnNonZero[jj][ii]   += 1;
						isColumnNonZero[jj-1][ii] += 1;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj+1][ii-1] += 1;
						isColumnNonZero[jj][ii-1]   += 1;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj+1][ii+1] += 1;
						isColumnNonZero[jj][ii+1]   += 1;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj+2][ii] += 1;
						isColumnNonZero[jj+1][ii] += 1;
					}
				}
			}
			int numNonZeroColumns = 0;
			//std::cout << "(" << i << "," << j << ")|";
			for(int m=-2; m<=2; m++)
			{
				for(int l=-2; l<=2; l++)
				{
					//std::cout << isColumnNonZero[jj+m][ii+l] << ",";
					if(isColumnNonZero[jj+m][ii+l]) numNonZeroColumns++;
					isColumnNonZero[jj+m][ii+l] = 0;
				}
				//std::cout << "|";
			}
			//std::cout << numNonZeroColumns << std::endl;
			num_nonzeros += numNonZeroColumns;
		}
	}
	//std::cout << "Total nonzeros: " << num_nonzeros << std::endl;

	cusp::coo_matrix<int, real, host_memory> CHost(nx*ny, nx*ny, num_nonzeros);

	real valuesInColumns[5][5];

	for(int m=-2; m<=2; m++)
	{
		for(int l=-2; l<=2; l++)
		{
			valuesInColumns[jj+m][ii+l]=0.0;
		}
	}

	int idx = 0;
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			if(j>0)
			{
				int I = (j-1)*nx+i+N_u;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj-1][ii] += 1;
					valuesInColumns[jj-1][ii] -= 1.0;
					isColumnNonZero[jj][ii]   += 1;
					valuesInColumns[jj][ii]   += 1.0;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					real xi = DirectForcingSolver<memoryType>::coeffs[I];
					if(diff < -1)
					{
						isColumnNonZero[jj-2][ii] += 1;
						valuesInColumns[jj-2][ii] -= xi;
						isColumnNonZero[jj-1][ii] += 1;
						valuesInColumns[jj-1][ii] += xi;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj-1][ii-1] += 1;
						valuesInColumns[jj-1][ii-1] -= xi;
						isColumnNonZero[jj][ii-1]   += 1;
						valuesInColumns[jj][ii-1]   += xi;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj-1][ii+1] += 1;
						valuesInColumns[jj-1][ii+1] -= xi;
						isColumnNonZero[jj][ii+1]   += 1;
						valuesInColumns[jj][ii+1]   += xi;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj][ii]   += 1;
						valuesInColumns[jj][ii]   -= xi;
						isColumnNonZero[jj+1][ii] += 1;
						valuesInColumns[jj+1][ii] += xi;
					}
				}
			}
			if(i>0)
			{
				int I = j*(nx-1)+(i-1);
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj][ii-1] += 1;
					valuesInColumns[jj][ii-1] -= 1.0;
					isColumnNonZero[jj][ii]   += 1;
					valuesInColumns[jj][ii]   += 1.0;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					real xi = DirectForcingSolver<memoryType>::coeffs[I];
					if(diff < -1)
					{
						isColumnNonZero[jj-1][ii-1] += 1;
						valuesInColumns[jj-1][ii-1] -= xi;
						isColumnNonZero[jj-1][ii]   += 1;
						valuesInColumns[jj-1][ii]   += xi;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj][ii-2] += 1;
						valuesInColumns[jj][ii-2] -= xi;
						isColumnNonZero[jj][ii-1] += 1;
						valuesInColumns[jj][ii-1] += xi;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj][ii]   += 1;
						valuesInColumns[jj][ii]   -= xi;
						isColumnNonZero[jj][ii+1] += 1;
						valuesInColumns[jj][ii+1] += xi;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj+1][ii-1] += 1;
						valuesInColumns[jj+1][ii-1] -= xi;
						isColumnNonZero[jj+1][ii]   += 1;
						valuesInColumns[jj+1][ii]   += xi;
					}
				}
			}
			if(i<nx-1)
			{
				int I = j*(nx-1)+i;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj][ii+1] += 1;
					valuesInColumns[jj][ii+1] -= 1.0;
					isColumnNonZero[jj][ii]   += 1;
					valuesInColumns[jj][ii]   += 1.0;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					real xi = DirectForcingSolver<memoryType>::coeffs[I];
					if(diff < -1)
					{
						isColumnNonZero[jj-1][ii+1] += 1;
						valuesInColumns[jj-1][ii+1] -= xi;
						isColumnNonZero[jj-1][ii]   += 1;
						valuesInColumns[jj-1][ii]   += xi;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj][ii]   += 1;
						valuesInColumns[jj][ii]   -= xi;
						isColumnNonZero[jj][ii-1] += 1;
						valuesInColumns[jj][ii-1] += xi;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj][ii+2] += 1;
						valuesInColumns[jj][ii+2] -= xi;
						isColumnNonZero[jj][ii+1] += 1;
						valuesInColumns[jj][ii+1] += xi;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj+1][ii+1] += 1;
						valuesInColumns[jj+1][ii+1] -= xi;
						isColumnNonZero[jj+1][ii]   += 1;
						valuesInColumns[jj+1][ii]   += xi;
					}
				}
			}
			if(j<ny-1)
			{
				int I = j*nx+i+N_u;
				if(DirectForcingSolver<memoryType>::tags[I]==-1)
				{
					isColumnNonZero[jj+1][ii] += 1;
					valuesInColumns[jj+1][ii] -= 1.0;
					isColumnNonZero[jj][ii]   += 1;
					valuesInColumns[jj][ii]   += 1.0;
				}
				else
				{
					int diff = DirectForcingSolver<memoryType>::tags[I]-I;
					real xi = DirectForcingSolver<memoryType>::coeffs[I];
					if(diff < -1)
					{
						isColumnNonZero[jj][ii]   += 1;
						valuesInColumns[jj][ii]   -= xi;
						isColumnNonZero[jj-1][ii] += 1;
						valuesInColumns[jj-1][ii] += xi;
					}
					else if(diff == -1)
					{
						isColumnNonZero[jj+1][ii-1] += 1;
						valuesInColumns[jj+1][ii-1] -= xi;
						isColumnNonZero[jj][ii-1]   += 1;
						valuesInColumns[jj][ii-1]   += xi;
					}
					else if(diff == 1)
					{
						isColumnNonZero[jj+1][ii+1] += 1;
						valuesInColumns[jj+1][ii+1] -= xi;
						isColumnNonZero[jj][ii+1]   += 1;
						valuesInColumns[jj][ii+1]   += xi;
					}
					else if(diff > 1)
					{
						isColumnNonZero[jj+2][ii] += 1;
						valuesInColumns[jj+2][ii] -= xi;
						isColumnNonZero[jj+1][ii] += 1;
						valuesInColumns[jj+1][ii] += xi;
					}
				}
			}
			int row = j*nx+i;
			for(int m=-2; m<=2; m++)
			{
				for(int l=-2; l<=2; l++)
				{
					if(isColumnNonZero[jj+m][ii+l])
					{
						CHost.row_indices[idx] = row;
						CHost.column_indices[idx] = row + m*nx + l;
						CHost.values[idx] = valuesInColumns[jj+m][ii+l];
						if(CHost.row_indices[idx]==(ny/2)*nx+nx/2 && CHost.row_indices[idx]==CHost.column_indices[idx])
						{
							CHost.values[idx]+=CHost.values[idx];
						}
						idx++;
					}
					isColumnNonZero[jj+m][ii+l] = 0;
					valuesInColumns[jj+m][ii+l] = 0.0;
				}
			}
		}
	}
	CHost.sort_by_row_and_column();
	CHost.values[0] += CHost.values[0];
	NavierStokesSolver<memoryType>::C = CHost;
	//cusp::io::write_matrix_market_file(NavierStokesSolver<memoryType>::C, "C-generateC.mtx");
	cusp::blas::scal(NavierStokesSolver<memoryType>::C.values, dt);
}
*/

template class DFImprovedSolver<host_memory>;
template class DFImprovedSolver<device_memory>;
