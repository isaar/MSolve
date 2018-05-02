using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Tests
{
    static class ReorderingTest
    {
        public static void TestAMDNonZeros(Model2D model)
        {
            var dofEnumerator1 = DOFEnumeratorInterleaved.Create(model);
            var assembler = new GlobalDOKAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, dofEnumerator1);
            using (var factor = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                Console.WriteLine("Interleaved ordering: Non zeros in factorization = " + factor.NumNonZeros);
            }
            (int[] permutation, ReorderingStatistics stats) = Kuu.Reorder(new OrderingAMD());
            Console.WriteLine("Interleaved ordering + AMD: Non zeros in factorization = " + stats.FactorizedNumNonZeros);

            var dofEnumerator2 = DOFEnumeratorSeparate.Create(model);
            (Kuu, Kuc) = assembler.BuildGlobalMatrix(model, dofEnumerator2);
            using (var factor = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                Console.WriteLine("Separate ordering: Non zeros in factorization = " + factor.NumNonZeros);
            }
            (permutation, stats) = Kuu.Reorder(new OrderingAMD());
            Console.WriteLine("Separate ordering + AMD: Non zeros in factorization = " + stats.FactorizedNumNonZeros);

        }
    }
}
