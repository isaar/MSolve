using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
	public interface IShellSectionMaterial:IFiniteElementMaterial
	{
		new IShellSectionMaterial Clone();
		double[] MembraneForces { get; }
		double[] Moments { get; }
		IMatrixView MembraneConstitutiveMatrix { get; }
		IMatrixView BendingConstitutiveMatrix { get; }
		IMatrixView CouplingConstitutiveMatrix { get; }
		void UpdateMaterial(double[] membraneStrains, double[] bendingStrains);
	}
}