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
    /// Tests for <see cref="SemidefiniteCholeskySkyline"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SemidefiniteCholeskySkylineTests
    {
        [Fact]
        private static void TestFactorizationOfBeam2DStiffness()
        {
            double tolerance = 1E-9;

            // Unconstrained beam
            var beam2DFreeK = Matrix.CreateFromArray(Beam2DElementMatrix.UnconstrainedStiffness);
            SemidefiniteCholeskySkyline factorFreeK =
                SkylineMatrix.CreateFromMatrix(beam2DFreeK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorFreeK.DependentColumns.Count == 3); // 2D problem: 3 rigid body modes, 0 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(beam2DFreeK, factorFreeK.NullSpaceBasis, tolerance);

            // Pinned beam
            var beam2DPinnedK = Matrix.CreateFromArray(Beam2DElementMatrix.PinnedStiffness);
            SemidefiniteCholeskySkyline factorPinnedK =
                SkylineMatrix.CreateFromMatrix(beam2DPinnedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorPinnedK.DependentColumns.Count == 1); // 2D problem: 3 rigid body modes, 2 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(beam2DPinnedK, factorPinnedK.NullSpaceBasis, tolerance);

            // Clamped beam
            var beam2DClampedK = Matrix.CreateFromArray(Beam2DElementMatrix.ClampedStiffness);
            SemidefiniteCholeskySkyline factorClampedK =
                SkylineMatrix.CreateFromMatrix(beam2DClampedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 2D problem: 3 rigid body modes, 3 are constrained
        }

        [Fact]
        private static void TestFactorizationOfBeam3DStiffness()
        {
            double tolerance = 1E-6;

            // Unconstrained beam
            var beam3DFreeK = Matrix.CreateFromArray(Beam3DElementMatrix.UnconstrainedStiffness);
            SemidefiniteCholeskySkyline factorFreeK = 
                SkylineMatrix.CreateFromMatrix(beam3DFreeK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorFreeK.DependentColumns.Count == 6); // 3D problem: 6 rigid body modes, 0 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(beam3DFreeK, factorFreeK.NullSpaceBasis, tolerance);

            // Pinned beam
            var beam3DPinnedK = Matrix.CreateFromArray(Beam3DElementMatrix.PinnedStiffness);
            SemidefiniteCholeskySkyline factorPinnedK =
                SkylineMatrix.CreateFromMatrix(beam3DPinnedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorPinnedK.DependentColumns.Count == 3); // 3D problem: 6 rigid body modes, 3 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(beam3DPinnedK, factorPinnedK.NullSpaceBasis, tolerance);

            // Clamped beam
            var beam3DClampedK = Matrix.CreateFromArray(Beam3DElementMatrix.ClampedStiffness);
            SemidefiniteCholeskySkyline factorClampedK =
                SkylineMatrix.CreateFromMatrix(beam3DClampedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 3D problem: 6 rigid body modes, 6 are constrained
        }

        [Fact]
        private static void TestFactorizationOfQuad4Stiffness()
        {
            double tolerance = 1E-7;

            // Unconstrained
            var quad4FreeK = Matrix.CreateFromArray(Quad4ElementMatrix.UnconstrainedStiffness);
            SemidefiniteCholeskySkyline factorFreeK =
                SkylineMatrix.CreateFromMatrix(quad4FreeK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorFreeK.DependentColumns.Count == 3); // 2D problem: 3 rigid body modes, 0 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(quad4FreeK, factorFreeK.NullSpaceBasis, tolerance);

            // Pinned
            var quad4PinnedK = Matrix.CreateFromArray(Quad4ElementMatrix.PinnedStiffness);
            SemidefiniteCholeskySkyline factorPinnedK =
                SkylineMatrix.CreateFromMatrix(quad4PinnedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorPinnedK.DependentColumns.Count == 1); // 2D problem: 3 rigid body modes, 2 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(quad4PinnedK, factorPinnedK.NullSpaceBasis, tolerance);

            // Clamped
            var quad4ClampedK = Matrix.CreateFromArray(Quad4ElementMatrix.ClampedStiffness);
            SemidefiniteCholeskySkyline factorClampedK =
                SkylineMatrix.CreateFromMatrix(quad4ClampedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 2D problem: 3 rigid body modes, 3 are constrained
        }

        [Fact]
        private static void TestFactorizationOfHexa8Stiffness()
        {
            double tolerance = 1E-2;
            double toleranceAbaqus = 1E-4;

            // Unconstrained
            var hexa8FreeK = Matrix.CreateFromArray(Hexa8ElementMatrix.UnconstrainedStiffness);
            SemidefiniteCholeskySkyline factorFreeK =
                SkylineMatrix.CreateFromMatrix(hexa8FreeK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorFreeK.DependentColumns.Count == 6); // 3D problem: 6 rigid body modes, 0 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(hexa8FreeK, factorFreeK.NullSpaceBasis, tolerance);

            // Pinned
            var hexa8PinnedK = Matrix.CreateFromArray(Hexa8ElementMatrix.PinnedStiffness);
            SemidefiniteCholeskySkyline factorPinnedK =
                SkylineMatrix.CreateFromMatrix(hexa8PinnedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorPinnedK.DependentColumns.Count == 3); // 3D problem: 6 rigid body modes, 3 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(hexa8PinnedK, factorPinnedK.NullSpaceBasis, tolerance);

            // Clamped
            var hexa8ClampedK = Matrix.CreateFromArray(Hexa8ElementMatrix.ClampedStiffness);
            SemidefiniteCholeskySkyline factorClampedK =
                SkylineMatrix.CreateFromMatrix(hexa8ClampedK).FactorSemidefiniteCholesky(true, tolerance);
            Assert.True(factorClampedK.DependentColumns.Count == 0); // 3D problem: 6 rigid body modes, 6 are constrained

            // Abaqus unconstrained
            var hexa8AbaqusK = Matrix.CreateFromArray(Hexa8ElementMatrixAbaqus.UnconstrainedStiffness);
            SemidefiniteCholeskySkyline factorAbaqusK =
                SkylineMatrix.CreateFromMatrix(hexa8AbaqusK).FactorSemidefiniteCholesky(true, toleranceAbaqus);
            Assert.True(factorAbaqusK.DependentColumns.Count == 6); // 3D problem: 6 rigid body modes, 0 are constrained
            SemidefiniteCholeskyFullTests.CheckNullSpace(hexa8AbaqusK, factorAbaqusK.NullSpaceBasis, toleranceAbaqus);
        }
    }
}
