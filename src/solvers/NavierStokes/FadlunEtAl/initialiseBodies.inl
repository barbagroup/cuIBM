template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialiseBodies()
{
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	B.initialise(db, *NavierStokesSolver<memoryType>::domInfo);
	std::cout << "Initialised bodies!" << std::endl;
}
