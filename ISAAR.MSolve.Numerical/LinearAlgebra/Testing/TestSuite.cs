using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestLibs;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing
{
    // TODO: need a symmetric, invertible, indefinite benchmark.
    public static class TestSuite
    {
        public static void TestAll()
        {
            SuiteSparse();
            //TestMarshaling();
            //TestBuilders();
            //TestInverse();
            //TestIndexing();
            //TestEquality();
            //TestFactorization();
            //TestMatrixOperations();
            //TestMatrixVectorMultiplication();
            //TestMatrixMatrixMultiplication();
            //TestSystemSolution();
            //TestTransposition();
            //TestVectorOperations();
            //TestWriting();
        }

        public static void TestMarshaling()
        {
            TestMKL.TestDgemv();
            TestMKL.TestDgetrf_Dgetrs();
        }

        public static void TestBuilders()
        {
            SparseRect.CheckBuilders();
        }

        public static void TestEquality()
        {
            //SparseRect.CheckEquals();
            SparsePositiveDefinite.CheckEquals();
        }

        public static void TestFactorization()
        {
            // Triangulations
            //SquareInvertible.CheckFactorization();
            //SquareSingular.CheckFactorization();
            //SquareSingular1Deficiency.CheckFactorization();
            //SymmPositiveDefinite.CheckFactorization();
            //SymmSingular.CheckFactorization();
            //SparsePositiveDefinite.CheckFactorization();

            // Orthogonalizations
            RectangularFullColRank.CheckFactorizationQR();
            RectangularFullColRank.CheckFactorizationLQ();
        }

        public static void TestIndexing()
        {
            //LowerInvertible.CheckIndexing();
            //LowerSingular.CheckIndexing();
            //UpperInvertible.CheckIndexing();
            //UpperSingular.CheckIndexing();
            //SymmPositiveDefinite.CheckIndexing();
            //SymmSingular.CheckIndexing();
            //SparseRect.CheckIndexing();
            SparsePositiveDefinite.CheckIndexing();
        }

        public static void TestInverse()
        {
            SquareInvertible.CheckInverse();
            SymmPositiveDefinite.CheckInverse();
        }

        public static void TestWriting()
        {
            DenseVectors.Print();
            Console.WriteLine();
            RectangularFullColRank.Print();
            Console.WriteLine();
            SquareInvertible.Print();
            Console.WriteLine();
            SquareSingular.Print();
            Console.WriteLine();
            SquareSingular1Deficiency.Print();
            Console.WriteLine();
            LowerInvertible.Print();
            Console.WriteLine();
            LowerSingular.Print();
            Console.WriteLine();
            UpperInvertible.Print();
            Console.WriteLine();
            UpperSingular.Print();
            Console.WriteLine();
            SymmPositiveDefinite.Print();
            Console.WriteLine();
            SymmSingular.Print();
            Console.WriteLine();
            SparseRect.Print();
            SparsePositiveDefinite.Print();
            Console.WriteLine();
        }

        public static void TestMatrixVectorMultiplication()
        {
            //Rectangular.CheckMatrixVectorMult();
            //SquareInvertible.CheckMatrixVectorMult();
            //SquareSingular.CheckMatrixVectorMult();
            //SquareSingular1Deficiency.CheckMatrixVectorMult();
            //LowerInvertible.CheckMatrixVectorMult();
            //LowerSingular.CheckMatrixVectorMult();
            //UpperInvertible.CheckMatrixVectorMult();
            //UpperSingular.CheckMatrixVectorMult();
            //SymmPositiveDefinite.CheckMatrixVectorMult();
            //SymmSingular.CheckMatrixVectorMult();
            //SparseRect.CheckMatrixVectorMult();
            SparsePositiveDefinite.CheckMatrixVectorMult();
        }

        public static void TestMatrixMatrixMultiplication()
        {
            SparseRect.CheckMatrixMatrixMult();
        }

        public static void TestSystemSolution()
        {
            // Linear systems
            //SquareInvertible.CheckSystemSolution();
            //SquareSingular.CheckSystemSolution();
            //SquareSingular1Deficiency.CheckSystemSolution();
            //LowerInvertible.CheckSystemSolution();
            //LowerSingular.CheckSystemSolution();
            //UpperInvertible.CheckSystemSolution();
            //UpperSingular.CheckSystemSolution();
            //SymmPositiveDefinite.CheckSystemSolution();
            //SymmSingular.CheckSystemSolution();
            //SparsePositiveDefinite.CheckSystemSolution();

            // Least squares systems
            RectangularFullColRank.CheckSolutionLeastSquares();
            RectangularFullColRank.CheckSolutionMinNorm();
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

        public static void SuiteSparse()
        {
            TestSuiteSparse.ExampleRawArrays();
            TestSuiteSparse.ExampleMatrixClasses();
            TestSuiteSparse.CheckRowAddition();
            //TestSuiteSparse.CheckRowDeletion();
        }
    }
}
