template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialiseBodies()
{
  parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	B.initialise(db, *NavierStokesSolver<memoryType>::domInfo);
}
