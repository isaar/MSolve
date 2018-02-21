using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing
{
    public static class TestSuite
    {
        public static void TestAll()
        {
            TestIndexing();
            //TestMatrixVectorMultiplication();
            //TestFactorization();
            //TestSystemSolution();
        }

        public static void TestFactorization()
        {
            SquareInvertible.CheckFactorization();
            SquareSingular.CheckFactorization();
            SquareSingular1Deficiency.CheckFactorization();
        }

        public static void TestIndexing()
        {
            LowerInvertible.CheckIndexing();
            LowerSingular.CheckIndexing();
            UpperInvertible.CheckIndexing();
            UpperSingular.CheckIndexing();
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
        }

        

    }
}
