using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData.FiniteElementMatrices
{
    /// <summary>
    /// Stiffness matrix of 2D Euler beam element with axial deformation (also called frame element). Various boundary 
    /// conditions are possible.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class Beam2DElementMatrix
    {
        internal static double[,] UnconstrainedStiffness
        {
            get
            {
                double E = 2.0E7, L = 2.0, b = 0.2, h = 0.4;
                double A = b * h;
                double I = b * h*h*h / 12.0;
                return new double[,] {
                    { E*A/L, 0, 0, -E*A/L, 0, 0 },
                    { 0, 12 * E * I / (L*L*L), 6 * E * I / (L*L), 0, -12 * E * I / (L*L*L), 6 * E * I / (L*L) },
                    { 0, 6 * E * I / (L*L), 4 * E * I / L, 0, -6 * E * I / (L*L), 2 * E * I / L },
                    { -E * A / L, 0, 0, E* A/ L, 0, 0 },
                    { 0, -12 * E * I / (L*L*L), -6 * E * I / (L*L), 0, 12 * E * I / (L*L*L), -6 * E * I / (L*L) },
                    { 0, 6 * E * I / (L*L), 2 * E * I / L, 0, -6 * E * I / (L*L), 4 * E * I / L }
                };
            }
        }

        internal static double[,] ClampedStiffness
        {
            get
            {
                var unconstrained = Matrix.CreateFromArray(UnconstrainedStiffness);
                int order = unconstrained.NumColumns;

                // constrain all dofs of the first node
                return unconstrained.GetSubmatrix(3, order, 3, order).CopyToArray2D(); 
            }
        }

        internal static double[,] PinnedStiffness
        {
            get
            {
                var unconstrained = Matrix.CreateFromArray(UnconstrainedStiffness);
                int order = unconstrained.NumColumns;

                // constrain translational dofs of the first node
                return unconstrained.GetSubmatrix(2, order, 2, order).CopyToArray2D(); 
            }
        }
    }
}
