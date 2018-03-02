using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing
{
    // TODO: need a symmetric, invertible, indefinite benchmark.
    public static class TestSuite
    {
        public static void TestAll()
        {
            //TestIndexing();
            TestMatrixVectorMultiplication();
            //TestFactorization();
            //TestSystemSolution();
            //TestVectorOperations();
            //TestPrinting();
        }

        public static void TestFactorization()
        {
            //SquareInvertible.CheckFactorization();
            //SquareSingular.CheckFactorization();
            //SquareSingular1Deficiency.CheckFactorization();
            SymmPositiveDefinite.CheckFactorization();
            SymmSingular.CheckFactorization();
        }

        public static void TestIndexing()
        {
            //LowerInvertible.CheckIndexing();
            //LowerSingular.CheckIndexing();
            //UpperInvertible.CheckIndexing();
            //UpperSingular.CheckIndexing();
            SymmPositiveDefinite.CheckIndexing();
            SymmSingular.CheckIndexing();
        }

        public static void TestPrinting()
        {
            DenseVectors.Print();
            SquareInvertible.Print();
        }

        public static void TestMatrixVectorMultiplication()
        {
            SquareInvertible.CheckMatrixVectorMult();
            SquareSingular.CheckMatrixVectorMult();
            SquareSingular1Deficiency.CheckMatrixVectorMult();
            LowerInvertible.CheckMatrixVectorMult();
            LowerSingular.CheckMatrixVectorMult();
            UpperInvertible.CheckMatrixVectorMult();
            UpperSingular.CheckMatrixVectorMult();
            SymmPositiveDefinite.CheckMatrixVectorMult();
            SymmSingular.CheckMatrixVectorMult();
        }

        public static void TestSystemSolution()
        {
            SquareInvertible.CheckSystemSolution();
            SquareSingular.CheckSystemSolution();
            SquareSingular1Deficiency.CheckSystemSolution();
            LowerInvertible.CheckSystemSolution();
            LowerSingular.CheckSystemSolution();
            UpperInvertible.CheckSystemSolution();
            UpperSingular.CheckSystemSolution();
            SymmPositiveDefinite.CheckSystemSolution();
            SymmSingular.CheckSystemSolution();
        }

        public static void TestVectorOperations()
        {
            DenseVectors.CheckAddition();
            DenseVectors.CheckAxpy();
            DenseVectors.CheckDotProduct();
            DenseVectors.CheckHadamardProduct();
            DenseVectors.CheckLinearCombination();
            DenseVectors.CheckNorm2();
            DenseVectors.CheckScaling();
            DenseVectors.CheckSubtraction();
        }

    }
}
