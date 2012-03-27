template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialiseBodies()
{
	B.initialise(*NavierStokesSolver<memoryType>::flowDesc, *NavierStokesSolver<memoryType>::domInfo);
}