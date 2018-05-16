using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Use a dedicated DiagonalMatrix class, instead of passing in double[] or Vector. It will also implement the inverse and 
// multiplication routines.
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning
{
    public class JacobiPreconditioner: IPreconditioner
    {
        public const double Tolerance = 1e-10;
        private readonly double[] inverseDiagonal;

        public JacobiPreconditioner(double[] diagonal, double tolerance = Tolerance)
        {
            Order = diagonal.Length;
            inverseDiagonal = new double[Order];
            for (int i = 0; i < Order; ++i)
            {
                double val = diagonal[i];
                if (Math.Abs(val) < tolerance) throw new SingularMatrixException($"Zero diagonal entry at index {i}");
                inverseDiagonal[i] = 1.0 / val;
            }
        }

        public int Order { get; }

        public Vector SolveLinearSystem(Vector rhs)
        {
            Preconditions.CheckSystemSolutionDimensions(Order, Order, rhs.Length);
            double[] solution = new double[Order];
            for (int i = 0; i < Order; ++i) solution[i] = inverseDiagonal[i] * rhs[i];
            return Vector.CreateFromArray(solution);
        }
    }
}
