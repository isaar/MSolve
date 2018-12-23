using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class CholeskySuiteSparseSolver: SolverBase
    {
        //private static readonly string directoryPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName;
        //private static readonly string matrixPath = directoryPath + "\\Resources\\matrixCSC.txt";
        //private static readonly string rhsPath = directoryPath + "\\Resources\\rhs.txt";
        //private static readonly string solutionPath = directoryPath + "\\Resources\\solution.txt";

        public CholeskySuiteSparseSolver(Model2D model) : base(model)
        { }

        public override void Solve()
        {
            DofOrderer = InterleavedDofOrderer.Create(model);
            //DofOrderer = DofOrdererSeparate.Create(model);
            var assembler = new GlobalDOKAssembler();
            (DokSymmetric Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            //if (!Kuu.IsSymmetric(double.Epsilon)) throw new AsymmetricMatrixException(
            //    "Stiffness matrix corresponding to free-free dofs is not symmetric");
            Vector rhs = CalcEffectiveRhs(Kuc);

            #region debug
            //// Matrix
            //(double[] values, int[] rowIndices, int[] colOffsets) = Kuu.BuildSymmetricCSCArrays(true);
            //var csc = CSCMatrix.CreateFromArrays(Kuu.NumRows, Kuu.NumColumns, values, rowIndices, colOffsets, false);
            //(new RawArraysWriter(csc)).WriteToMultipleFiles(matrixPath, true);

            //// Rhs
            //FullVectorWriter.NumericFormat = new GeneralNumericFormat();
            //(new FullVectorWriter(rhs, true)).WriteToFile(rhsPath);

            //// Solution
            //Vector expectedSolution = SolveWithSkyline();
            //(new FullVectorWriter(expectedSolution, true)).WriteToFile(solutionPath);
            #endregion

            using (var factor = CholeskySuiteSparse.Factorize(Kuu.BuildSymmetricCscMatrix(true), true))
            {
                Solution = factor.SolveLinearSystem(rhs);
            }
        }

        private Vector SolveWithSkyline()
        {
            var assembler = new GlobalSkylineAssembler();
            (SkylineMatrix Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);
            return Kuu.FactorCholesky(true).SolveLinearSystem(rhs);
        }
    }
}
