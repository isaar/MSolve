using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Materials.Interfaces
{
	public interface IShellSectionMaterial:IFiniteElementMaterial
	{
		new IShellSectionMaterial Clone();
		double[] MembraneForces { get; }
		double[] Moments { get; }
		IMatrix2D MembraneConstitutiveMatrix { get; }
		IMatrix2D BendingConstitutiveMatrix { get; }
		IMatrix2D CouplingConstitutiveMatrix { get; }
		void UpdateMaterial(double[] membraneStrains, double[] bendingStrains);
	}
}