using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;

//TODO: probably needs a builder
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1
{
    /// <summary>
    /// Uses the standard Preconditioned Conjugate Gradient to solve the interface linear system after projecting it to a new 
    /// space. The particular solution of the lagrange multipliers is calculated separately and at the end is added to the PCG
    /// solution.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Feti1ProjectedInterfaceProblemSolver : IFeti1InterfaceProblemSolver
    {
        public enum ProjectionSide
        {
            None, Left, Both
        }

        private readonly double pcgConvergenceTolerance;
        private readonly double pcgMaxIterationsOverSize;
        private readonly ProjectionSide matrixProjectionSide;
        private readonly ProjectionSide preconditionerProjectionSide;

        public Feti1ProjectedInterfaceProblemSolver(double pcgConvergenceTolerance, double pcgMaxIterationsOverSize,
            ProjectionSide matrixProjectionSide, ProjectionSide preconditionerProjectionSide)
        {
            this.pcgConvergenceTolerance = pcgConvergenceTolerance;
            this.pcgMaxIterationsOverSize = pcgMaxIterationsOverSize;
            this.matrixProjectionSide = matrixProjectionSide;
            this.preconditionerProjectionSide = preconditionerProjectionSide;
        }

        public Vector CalcLagrangeMultipliers(Feti1FlexibilityMatrix flexibility, IFetiPreconditioner preconditioner, 
            Feti1Projection projection, Vector disconnectedDisplacements, Vector rigidBodyModesWork, double globalForcesNorm,
            FetiLogger logger)
        {
            int systemOrder = flexibility.Order;

            PcgMatrix pcgMatrix = DefinePcgMatrix(flexibility, projection);
            PcgPreconditioner pcgPreconditioner = DefinePcgPreconditioner(preconditioner, projection);

            // λ0 = Q * G * inv(G^T * Q * G) * e
            Vector lagrangesParticular = projection.CalcParticularLagrangeMultipliers(rigidBodyModesWork);

            // Calculate rhs of the projected interface system: rhs = P^T * (d - F * λ0)
            var r0 = flexibility.Multiply(lagrangesParticular);
            r0.LinearCombinationIntoThis(-1.0, disconnectedDisplacements, 1.0);
            var pcgRhs = Vector.CreateZero(systemOrder);
            projection.ProjectVector(r0, pcgRhs, true);

            // Solve the interface problem using PCG algorithm
            var pcgBuilder = new PcgAlgorithm.Builder();
            pcgBuilder.MaxIterationsProvider = new PercentageMaxIterationsProvider(pcgMaxIterationsOverSize);
            pcgBuilder.ResidualTolerance = pcgConvergenceTolerance;
            //pcgBuilder.ResidualNormCalculator = calcExactResidualNorm(); //TODO: implement this for regular PCG too. WARNING. λ = λ0 + P^T * λbar is needed.
            PcgAlgorithm pcg = pcgBuilder.Build();
            var lagrangesBar = Vector.CreateZero(systemOrder);
            CGStatistics stats = pcg.Solve(pcgMatrix, pcgPreconditioner, pcgRhs, lagrangesBar, true,
                () => Vector.CreateZero(systemOrder));

            // Log statistics about PCG execution
            if (!stats.HasConverged)
            {
                throw new IterativeSolverNotConvergedException(Feti1Solver.name + " did not converge to a solution. PCG"
                    + $" algorithm run for {stats.NumIterationsRequired} iterations and the residual norm ratio was"
                    + $" {stats.NormRatio}");
            }
            logger.PcgIterations = stats.NumIterationsRequired;
            logger.PcgResidualNormRatio = stats.NormRatio;

            // Calculate the lagrange multipliers from the separation formula: λ = λ0 + P * λbar
            var lagranges = Vector.CreateZero(systemOrder);
            projection.ProjectVector(lagrangesBar, lagranges, false);
            lagranges.AddIntoThis(lagrangesParticular);
            return lagranges;
        }

        private PcgMatrix DefinePcgMatrix(Feti1FlexibilityMatrix flexibility, Feti1Projection projection)
        {
            //TODO: projection.ProjectVector() overwrites the rhs vector passed in. However that is not clear by the method.
            //      I should either always return the resulting vector (negating any optimizations) or set the convention that
            //      vectors passed in as results will always be overwritten. The latter forces the method that accepts them to 
            //      clear them (if the method needs to), which is not needed if the resulting vector has just be initialized to 0.

            if (matrixProjectionSide == ProjectionSide.Left)
            {
                // A * x = (P^T * F) * x = P^T * (F * x)
                return new PcgMatrix(flexibility.Order, (x, y) => projection.ProjectVector(flexibility.Multiply(x), y, true));
            }
            else if (matrixProjectionSide == ProjectionSide.Both)
            {
                // A * x = (P^T * F * P) * x = P^T * (F * (P * x))
                return new PcgMatrix(flexibility.Order, (x, y) =>
                {
                    var temp = Vector.CreateZero(flexibility.Order);
                    projection.ProjectVector(x, temp, false);
                    projection.ProjectVector(flexibility.Multiply(temp), y, true);
                });
            }
            else throw new ArgumentException("The FETI-1 flexibility matrix can be multiplied from the left or both sides.");
        }

        private PcgPreconditioner DefinePcgPreconditioner(IFetiPreconditioner preconditioner, Feti1Projection projection)
        {
            if (matrixProjectionSide == ProjectionSide.None)
            {
                // x = prec(A) * y => x = prec(F) * y
                return new PcgPreconditioner((x, y) => preconditioner.SolveLinearSystem(y, x));
            }
            else if (matrixProjectionSide == ProjectionSide.Left)
            {
                // x = prec(A) * y => x = (P * prec(F)) * y => x = P * (prec(F) * y)
                return new PcgPreconditioner((y, x) =>
                {
                    int order = y.Length;
                    var temp = Vector.CreateZero(order);
                    preconditioner.SolveLinearSystem(y, temp);
                    projection.ProjectVector(temp, x, false);
                });
            }
            else if (matrixProjectionSide == ProjectionSide.Both)
            {
                // x = prec(A) * y => x = (P * prec(F) * P^T) * y => x = P * (prec(F) * (P^T * y))
                return new PcgPreconditioner((y, x) =>
                {
                    int order = y.Length;
                    var temp1 = Vector.CreateZero(order);
                    projection.ProjectVector(y, temp1, true);
                    var temp2 = Vector.CreateZero(order);
                    preconditioner.SolveLinearSystem(temp1, temp2);
                    projection.ProjectVector(temp2, x, false);
                });
            }
            else throw new ArgumentException(
                "The FETI-1 preconditioner can be multiplied from the left side, both sides or not at all.");
        }

        //TODO: Perhaps PCG, PCPG could accept delegates as well as matrix/preconditioner classes. Then client code needs to 
        //      define anonymous methods instead of classes.
        private class PcgMatrix : ILinearTransformation
        {
            private readonly Action<Vector, Vector> multiply;

            internal PcgMatrix(int order, Action<Vector, Vector> multiply)
            {
                this.NumColumns = order;
                this.multiply = multiply;
            }

            public int NumColumns { get; }
            public int NumRows => NumColumns;

            public void Multiply(IVectorView lhsVector, IVector rhsVector)
            {
                //TODO: remove casts. I think PCG, LinearTransformation and preconditioners should be generic, bounded by 
                //      IVectorView and IVector
                var lhs = (Vector)lhsVector;
                var rhs = (Vector)rhsVector;
                multiply(lhs, rhs);
            }
        }

        private class PcgPreconditioner : IPreconditioner
        {
            private readonly Action<Vector, Vector> precondition;

            internal PcgPreconditioner(Action<Vector, Vector> precondition)
            {
                this.precondition = precondition;
            }

            public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
            {
                //TODO: remove casts. I think PCG, LinearTransformation and preconditioners should be generic, bounded by 
                //      IVectorView and IVector
                var lhs = (Vector)lhsVector;
                var rhs = (Vector)rhsVector;
                precondition(rhs, lhs);
            }
        }
    }
}
