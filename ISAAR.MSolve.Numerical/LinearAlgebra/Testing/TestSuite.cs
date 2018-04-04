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
            //TestSuiteSparse.ExampleRawArrays();
            //TestSuiteSparse.ExampleMatrixClasses();
            //TestIndexing();
            //CheckEquality();
            //TestMatrixVectorMultiplication();
            TestTransposition();
            //TestFactorization();
            //TestSystemSolution();
            //TestVectorOperations();
            //TestMatrixOperations();
            //TestPrinting();
        }

        public static void CheckEquality()
        {
            SparseRect.CheckEquals();
        }

        public static void TestFactorization()
        {
            SquareInvertible.CheckFactorization();
            SquareSingular.CheckFactorization();
            SquareSingular1Deficiency.CheckFactorization();
            //SymmPositiveDefinite.CheckFactorization();
            //SymmSingular.CheckFactorization();
        }

        public static void TestIndexing()
        {
            LowerInvertible.CheckIndexing();
            LowerSingular.CheckIndexing();
            //UpperInvertible.CheckIndexing();
            //UpperSingular.CheckIndexing();
            //SymmPositiveDefinite.CheckIndexing();
            //SymmSingular.CheckIndexing();
            //SparseRect.CheckIndexing();
        }

        public static void TestPrinting()
        {
            DenseVectors.Print();
            SquareInvertible.Print();
        }

        public static void TestMatrixVectorMultiplication()
        {
            //Rectangular.CheckMatrixVectorMult();
            //SquareInvertible.CheckMatrixVectorMult();
            //SquareSingular.CheckMatrixVectorMult();
            //SquareSingular1Deficiency.CheckMatrixVectorMult();
            LowerInvertible.CheckMatrixVectorMult();
            LowerSingular.CheckMatrixVectorMult();
            //UpperInvertible.CheckMatrixVectorMult();
            //UpperSingular.CheckMatrixVectorMult();
            //SymmPositiveDefinite.CheckMatrixVectorMult();
            //SymmSingular.CheckMatrixVectorMult();
            //SparseRect.CheckMatrixVectorMult();
        }

        public static void TestSystemSolution()
        {
            //SquareInvertible.CheckSystemSolution();
            //SquareSingular.CheckSystemSolution();
            //SquareSingular1Deficiency.CheckSystemSolution();
            LowerInvertible.CheckSystemSolution();
            LowerSingular.CheckSystemSolution();
            //UpperInvertible.CheckSystemSolution();
            //UpperSingular.CheckSystemSolution();
            //SymmPositiveDefinite.CheckSystemSolution();
            //SymmSingular.CheckSystemSolution();
        }

        public static void TestTransposition()
        {
            DenseMatrices.CheckTransposition();
            LowerInvertible.CheckTransposition();
            LowerSingular.CheckTransposition();
            UpperInvertible.CheckTransposition();
            UpperSingular.CheckTransposition();
            SymmPositiveDefinite.CheckTransposition();
            SymmSingular.CheckTransposition();
            SparseRect.CheckTransposition();
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

        public static void TestMatrixOperations()
        {
            DenseMatrices.CheckScaling();
            DenseMatrices.CheckAddition();
            DenseMatrices.CheckSubtraction();
            DenseMatrices.CheckLinearCombination();
            DenseMatrices.CheckTransposition();
            DenseMatrices.CheckMatrixMultiplication();
        }

    }
}
