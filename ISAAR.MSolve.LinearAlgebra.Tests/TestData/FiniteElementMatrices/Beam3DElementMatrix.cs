using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData.FiniteElementMatrices
{
    /// <summary>
    /// Stiffness matrix of 3D Euler beam element with axial deformation (also called frame element). Various boundary 
    /// conditions are possible.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class Beam3DElementMatrix
    {
        internal static double[,] UnconstrainedStiffness
        {
            get
            {
                double E = 21000.0, poisson = 0.3, G = E / (2 + 2 * poisson);
                double L = 2.0;
                double A = 91.04;
                double Iy = 2843.0;
                double Iz = 8091.0;
                double J = 76.57;

                double L2 = L * L;
                double L3 = L * L * L;

                return new double[,] {
                    { E*A/L, 0, 0, 0, 0, 0, -E*A/L, 0, 0, 0, 0, 0 },
                    { 0, 12*E*Iz/L3, 0, 0, 0, 6*E*Iz/L2, 0, -12*E*Iz/L3, 0, 0, 0, 6*E*Iz/L2},
                    { 0, 0, 12*E*Iy/L3, 0, -6*E*Iy/L2, 0, 0, 0, -12*E*Iy/L3, 0, -6*E*Iy/L2, 0},
                    { 0, 0, 0, G*J/L, 0, 0, 0, 0, 0, -G*J/L, 0, 0 },
                    { 0, 0, -6*E*Iy/L2, 0, 4*E*Iy/L, 0, 0, 0, 6*E*Iy/L2, 0, 2*E*Iy/L, 0},
                    { 0, 6*E*Iz/L2, 0, 0, 0, 4*E*Iz/L, 0, -6*E*Iz/L2, 0, 0, 0, 2*E*Iz/L},
                    { -E*A/L, 0, 0, 0, 0, 0, E*A/L, 0, 0, 0, 0, 0 },
                    { 0, -12*E*Iz/L3, 0, 0, 0, -6*E*Iz/L2, 0, 12*E*Iz/L3, 0, 0, 0, -6*E*Iz/L2},
                    { 0, 0, -12*E*Iy/L3, 0, 6*E*Iy/L2, 0, 0, 0, 12*E*Iy/L3, 0, 6*E*Iy/L2, 0},
                    { 0, 0, 0, -G*J/L, 0, 0, 0, 0, 0, G*J/L, 0, 0 },
                    { 0, 0, -6*E*Iy/L2, 0, 2*E*Iy/L, 0, 0, 0, 6*E*Iy/L2, 0, 4*E*Iy/L, 0},
                    { 0, 6*E*Iz/L2, 0, 0, 0, 2*E*Iz/L, 0, -6*E*Iz/L2, 0, 0, 0, 4*E*Iz/L}
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
                return unconstrained.GetSubmatrix(6, order, 6, order).CopyToArray2D();
            }
        }

        internal static double[,] PinnedStiffness
        {
            get
            {
                var unconstrained = Matrix.CreateFromArray(UnconstrainedStiffness);
                int order = unconstrained.NumColumns;

                // constrain translational dofs of the first node
                return unconstrained.GetSubmatrix(3, order, 3, order).CopyToArray2D();
            }
        }
    }
}
