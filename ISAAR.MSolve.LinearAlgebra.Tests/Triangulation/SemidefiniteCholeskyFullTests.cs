using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData.FiniteElementMatrices;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

//TODO: Also check the invertible part of the factorized data 
namespace ISAAR.MSolve.LinearAlgebra.Tests.Triangulation
{
    /// <summary>
    /// Tests for <see cref="SemidefiniteCholeskyFull"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SemidefiniteCholeskyFullTests
    {
        internal static void CheckNullSpace(IMatrixView unfactorizedMatrix, IReadOnlyList<double[]> nullSpaceBasis, 
            double tolerance)
        {
            var comparer = new MatrixComparer(tolerance);

            // Check that each vector belongs to the nullspace
            int order = unfactorizedMatrix.NumColumns;
            //var zeroVector = Vector.CreateZero(order);
            int nullity = nullSpaceBasis.Count;
            var nullSpaceMatrix = Matrix.CreateZero(order, nullity);

            for (int j = 0; j < nullity; ++j)
            {
                var x = Vector.CreateFromArray(nullSpaceBasis[j]);
                nullSpaceMatrix.SetSubcolumn(j, x);

                // Check that each vector belongs to the nullspace
                IVector Ax = unfactorizedMatrix.Multiply(x);
                double normAx = Ax.Norm2() / Ax.Length;
                comparer.AssertEqual(0.0, normAx);
            }

            // Check that the vectors are independent. 
            // TODO: perhaps this should be included in the LinearAlgebra project, not just for tests.
            (Matrix rref, List<int> independentCols) = nullSpaceMatrix.ReducedRowEchelonForm(tolerance);
            for (int j = 0; j < nullity; ++j) Assert.Contains(j, independentCols);
        }

        [Fact]
        private static void TestFactorizationOfBeam2DStiffness()
        {
            double tolerance = 1E-9;

            // Unconstrained beam
            var beam2DFreeK = Matrix.CreateFromArray(Beam2DElementMatrix.UnconstrainedStiffness);
            var factorFreeK = SemidefiniteCholeskyFull.Factorize(beam2DFreeK.Copy(), tolerance);
            Assert.True(factorFreeK.DependentColumns.Count == 3); // 2D problem: 3 rigid body modes, 0 are constrained
            CheckNullSpace(beam2DFreeK, factorFreeK.NullSpaceBasis, tolerance);

            // Pinned beam
            var beam2DPinnedK = Matrix.CreateFromArray(Beam2DElementMatrix.PinnedStiffness);
            var factorPinnedK = SemidefiniteCholeskyFull.Factorize(beam2DPinnedK.Copy(), tolerance);
            Assert.True(factorPinnedK.DependentColumns.Count == 1); // 2D problem: 3 rigid body modes, 2 are constrained
            CheckNullSpace(beam2DPinnedK, factorPinnedK.NullSpaceBasis, tolerance);

            // Clamped beam
            var beam2DClampedK = Matrix.CreateFromArray(Beam2DElementMatrix.ClampedStiffness);
            var factorClampedK = SemidefiniteCholeskyFull.Factorize(beam2DClampedK.Copy(), tolerance);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 2D problem: 3 rigid body modes, 3 are constrained
        }

        [Fact]
        private static void TestFactorizationOfBeam3DStiffness()
        {
            double tolerance = 1E-6;

            // Unconstrained beam
            var beam3DFreeK = Matrix.CreateFromArray(Beam3DElementMatrix.UnconstrainedStiffness);
            var factorFreeK = SemidefiniteCholeskyFull.Factorize(beam3DFreeK.Copy(), tolerance);
            Assert.True(factorFreeK.DependentColumns.Count == 6); // 3D problem: 6 rigid body modes, 0 are constrained
            CheckNullSpace(beam3DFreeK, factorFreeK.NullSpaceBasis, tolerance);

            // Pinned beam
            var beam3DPinnedK = Matrix.CreateFromArray(Beam3DElementMatrix.PinnedStiffness);
            var factorPinnedK = SemidefiniteCholeskyFull.Factorize(beam3DPinnedK.Copy(), tolerance);
            Assert.True(factorPinnedK.DependentColumns.Count == 3); // 3D problem: 6 rigid body modes, 3 are constrained
            CheckNullSpace(beam3DPinnedK, factorPinnedK.NullSpaceBasis, tolerance);

            // Clamped beam
            var beam3DClampedK = Matrix.CreateFromArray(Beam3DElementMatrix.ClampedStiffness);
            var factorClampedK = SemidefiniteCholeskyFull.Factorize(beam3DClampedK.Copy(), tolerance);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 3D problem: 6 rigid body modes, 6 are constrained
        }

        [Fact]
        private static void TestFactorizationOfQuad4Stiffness()
        {
            double tolerance = 1E-7;

            // Unconstrained
            var quad4FreeK = Matrix.CreateFromArray(Quad4ElementMatrix.UnconstrainedStiffness);
            var factorFreeK = SemidefiniteCholeskyFull.Factorize(quad4FreeK.Copy(), tolerance);
            Assert.True(factorFreeK.DependentColumns.Count == 3); // 2D problem: 3 rigid body modes, 0 are constrained
            CheckNullSpace(quad4FreeK, factorFreeK.NullSpaceBasis, tolerance);

            // Pinned
            var quad4PinnedK = Matrix.CreateFromArray(Quad4ElementMatrix.PinnedStiffness);
            var factorPinnedK = SemidefiniteCholeskyFull.Factorize(quad4PinnedK.Copy(), tolerance);
            Assert.True(factorPinnedK.DependentColumns.Count == 1); // 2D problem: 3 rigid body modes, 2 are constrained
            CheckNullSpace(quad4PinnedK, factorPinnedK.NullSpaceBasis, tolerance);

            // Clamped
            var quad4ClampedK = Matrix.CreateFromArray(Quad4ElementMatrix.ClampedStiffness);
            var factorClampedK = SemidefiniteCholeskyFull.Factorize(quad4ClampedK.Copy(), tolerance);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 2D problem: 3 rigid body modes, 3 are constrained
        }

        [Fact]
        private static void TestFactorizationOfHexa8Stiffness()
        {
            double toleranceMSolve = 1E-2;
            double toleranceAbaqus = 1E-4;

            // Unconstrained
            var hexa8FreeK = Matrix.CreateFromArray(Hexa8ElementMatrix.UnconstrainedStiffness);
            var factorFreeK = SemidefiniteCholeskyFull.Factorize(hexa8FreeK.Copy(), toleranceMSolve);
            Assert.True(factorFreeK.DependentColumns.Count == 6); // 3D problem: 6 rigid body modes, 0 are constrained
            CheckNullSpace(hexa8FreeK, factorFreeK.NullSpaceBasis, toleranceMSolve);

            // Pinned
            var hexa8PinnedK = Matrix.CreateFromArray(Hexa8ElementMatrix.PinnedStiffness);
            var factorPinnedK = SemidefiniteCholeskyFull.Factorize(hexa8PinnedK.Copy(), toleranceMSolve);
            Assert.True(factorPinnedK.DependentColumns.Count == 3); // 3D problem: 6 rigid body modes, 3 are constrained
            CheckNullSpace(hexa8PinnedK, factorPinnedK.NullSpaceBasis, toleranceMSolve);

            // Clamped
            var hexa8ClampedK = Matrix.CreateFromArray(Hexa8ElementMatrix.ClampedStiffness);
            var factorClampedK = SemidefiniteCholeskyFull.Factorize(hexa8ClampedK.Copy(), toleranceMSolve);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 3D problem: 6 rigid body modes, 6 are constrained

            // Abaqus unconstrained
            var hexa8AbaqusK = Matrix.CreateFromArray(Hexa8ElementMatrixAbaqus.UnconstrainedStiffness);
            var factorAbaqusK = SemidefiniteCholeskyFull.Factorize(hexa8AbaqusK.Copy(), toleranceAbaqus);
            Assert.True(factorAbaqusK.DependentColumns.Count == 6); // 3D problem: 6 rigid body modes, 0 are constrained
            CheckNullSpace(hexa8AbaqusK, factorAbaqusK.NullSpaceBasis, toleranceAbaqus);
        }
    }
}
