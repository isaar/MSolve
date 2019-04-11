using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    public abstract class GmresAlgorithmBase
    {
        //public GmresStatistics Solve(IMatrixView matrix, IVector rhs, IVector solution,
        //    IPreconditioner preconditioner,int restartIterations, int maximumIterations, double tolerance)
        //{
        //    var iteration = 0;
        //    var flag = 0;

        //    var rhsNorm = rhs.Norm2();
        //    if (Math.Abs(rhsNorm) < 1e-9)
        //        rhsNorm = 1;

        //    IVector residualVector = Vector.CreateZero(rhs.Length);
        //    preconditioner.SolveLinearSystem(rhs.Subtract(matrix.Multiply(solution)),residualVector);

        //    var error = residualVector.Norm2() / rhsNorm;

        //    if (error < tolerance)
        //        return null;

            
        //}
    }
    
}
