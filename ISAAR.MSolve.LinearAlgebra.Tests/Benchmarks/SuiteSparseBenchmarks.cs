using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Benchmarks
{
    public class SuiteSparseBenchmarks
    {
        public static void MemoryConsumptionDebugging()
        {
            int order = 100000;
            int bandwidth = 200;
            for (int rep = 0; rep < 10; ++rep)
            {
                var dok = DokSymmetric.CreateEmpty(order);
                for (int i = 0; i < order; ++i)
                {
                    dok[i, i] = 10.0;
                    if (i >= bandwidth) dok[i - bandwidth, i] = 1.0;
                    else dok[0, i] = 1.0;
                }
                dok[0, 0] = 10.0;

                SymmetricCscMatrix matrix = dok.BuildSymmetricCscMatrix(true);
                var rhs = Vector.CreateWithValue(order, 2.0);
                var solution = Vector.CreateZero(order);

                using (var factorization = CholeskySuiteSparse.Factorize(matrix, true))
                {
                    factorization.SolveLinearSystem(rhs, solution);
                }
            }
        }
    }
}
