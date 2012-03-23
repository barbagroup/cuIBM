template <typename Matrix, typename Vector>
void TairaColoniusSolver<Matrix, Vector>::initialiseBodies()
{
	B.initialise(*NavierStokesSolver<Matrix, Vector>::flowDesc, *NavierStokesSolver<Matrix, Vector>::domInfo);
}