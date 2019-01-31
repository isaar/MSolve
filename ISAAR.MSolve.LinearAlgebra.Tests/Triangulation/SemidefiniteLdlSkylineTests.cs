using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Input;
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
    /// Tests for <see cref="SemidefiniteLdlSkyline"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class SemidefiniteLdlSkylineTests
    {
        private static void CheckFactorization(SkylineMatrix originalMatrix, int expectedNullity, double tolerance)
        {
            var factorization = originalMatrix.FactorSemidefiniteLdl(false, tolerance);
            Assert.True(factorization.DependentColumns.Count == expectedNullity,
                $"Nullity = {factorization.DependentColumns.Count} != {expectedNullity} (expected)");
            SemidefiniteCholeskyFullTests.CheckNullSpace(originalMatrix, factorization.NullSpaceBasis, tolerance);
        }

        [Fact]
        private static void TestFactorizationOfBeam2DStiffness()
        {
            double tolerance = 1E-11;

            // Unconstrained beam. 2D problem: 3 rigid body modes, 0 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Beam2DElementMatrix.UnconstrainedStiffness), 3, tolerance);

            // Pinned beam. 2D problem: 3 rigid body modes, 2 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Beam2DElementMatrix.PinnedStiffness), 1, tolerance);

            // Clamped beam. 2D problem: 3 rigid body modes, 3 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Beam2DElementMatrix.ClampedStiffness), 0, tolerance);
        }

        [Fact]
        private static void TestFactorizationOfBeam3DStiffness()
        {
            double tolerance = double.Epsilon;

            // Unconstrained beam. 3D problem: 6 rigid body modes, 0 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Beam3DElementMatrix.UnconstrainedStiffness), 6, tolerance);

            // Pinned beam. 3D problem: 6 rigid body modes, 3 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Beam3DElementMatrix.PinnedStiffness), 3, tolerance);

            // Clamped beam. 3D problem: 6 rigid body modes, 6 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Beam3DElementMatrix.ClampedStiffness), 0, tolerance);
        }

        [Fact]
        private static void TestFactorizationOfQuad4Stiffness()
        {
            double tolerance = 1E-7;

            // Unconstrained beam. 2D problem: 3 rigid body modes, 0 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Quad4ElementMatrix.UnconstrainedStiffness), 3, tolerance);

            // Pinned beam. 2D problem: 3 rigid body modes, 2 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Quad4ElementMatrix.PinnedStiffness), 1, tolerance);

            // Clamped beam. 2D problem: 3 rigid body modes, 3 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Quad4ElementMatrix.ClampedStiffness), 0, tolerance);
        }

        [Fact]
        private static void TestFactorizationOfHexa8Stiffness()
        {
            double toleranceMSolve = 1E-2;
            double toleranceAbaqus = 1E-4;

            // Unconstrained beam. 3D problem: 6 rigid body modes, 0 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Hexa8ElementMatrix.UnconstrainedStiffness), 6, toleranceMSolve);

            // Pinned beam. 3D problem: 6 rigid body modes, 3 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Hexa8ElementMatrix.PinnedStiffness), 3, toleranceMSolve);

            // Clamped beam. 3D problem: 6 rigid body modes, 6 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Hexa8ElementMatrix.ClampedStiffness), 0, toleranceMSolve);

            // Unconstrained beam. 3D problem: 6 rigid body modes, 0 are constrained
            CheckFactorization(SkylineMatrix.CreateFromArray(Hexa8ElementMatrixAbaqus.UnconstrainedStiffness), 6, toleranceAbaqus);
        }

        [Fact]
        private static void TestFactorizationOf2DStructureStiffness() // 2D problem: 3 rigid body modes, 0 are constrained
        {
            string resourcesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName + "\\Resources";
            string valuesPath = resourcesPath + "\\Quad4_20x20_stiffness_values.txt";
            string diagOffsetsPath = resourcesPath + "\\Quad4_20x20_stiffness_diagonal_offsets.txt";
            var reader = new RawArraysReader();
            SkylineMatrix unconstrainedK = reader.ReadSkylineMatrixFromSeparateFiles(valuesPath, diagOffsetsPath, true);

            double tolerance = 1E-2;
            CheckFactorization(unconstrainedK, 3, tolerance);
            CheckFactorization(unconstrainedK.Scale(1E-6), 3, 1E-8); //TODO: Normalization is required when comparing the pivot with the tolerance.
        }

        [Fact]
        private static void TestFactorizationOf3DStructureStiffness() // 3D problem: 6 rigid body modes, 0 are constrained
        {
            string resourcesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName + "\\Resources";
            string valuesPath = resourcesPath + "\\Hexa8_10x10x10_stiffness_values.txt";
            string diagOffsetsPath = resourcesPath + "\\Hexa8_10x10x10_stiffness_diagonal_offsets.txt";
            var reader = new RawArraysReader();
            SkylineMatrix unconstrainedK = reader.ReadSkylineMatrixFromSeparateFiles(valuesPath, diagOffsetsPath, true);

            double tolerance = 1E-3;
            CheckFactorization(unconstrainedK, 6, tolerance);
            CheckFactorization(unconstrainedK.Scale(1E-5), 6, 1E-8); //TODO: Normalization is required when comparing the pivot with the tolerance.
        }
    }
}
