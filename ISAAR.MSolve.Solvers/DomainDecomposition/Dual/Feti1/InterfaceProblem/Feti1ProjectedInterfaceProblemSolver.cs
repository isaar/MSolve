using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

//TODO: probably needs a builder
//TODO: all those enums are not objected oriented. Create strategies.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.InterfaceProblem
{
    /// <summary>
    /// Uses the standard Preconditioned Conjugate Gradient to solve the interface linear system after projecting it to a new 
    /// space. The particular solution of the lagrange multipliers is calculated separately and at the end is added to the PCG
    /// solution.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Feti1ProjectedInterfaceProblemSolver : IFeti1InterfaceProblemSolver
    {
        public enum LagrangeMultiplierSeparation
        {
            /// <summary>
            /// λ = λ0 + λbar
            /// </summary>
            Simple,

            /// <summary>
            /// λ = λ0 + P * λbar
            /// </summary>
            WithProjection
        }

        public enum ProjectionSide
        {
            None, Left, Both
        }
        
        private readonly LagrangeMultiplierSeparation lagrangeSeparation;
        private readonly IMaxIterationsProvider maxIterationsProvider;
        private readonly IFetiPcgConvergenceFactory pcgConvergenceStrategyFactory;
        private readonly double pcgConvergenceTolerance;
        private readonly ProjectionSide projectionSideMatrix;
        private readonly ProjectionSide projectionSidePreconditioner;

        private Feti1ProjectedInterfaceProblemSolver(IMaxIterationsProvider maxIterationsProvider, 
            double pcgConvergenceTolerance, IFetiPcgConvergenceFactory pcgConvergenceStrategyFactory,
            ProjectionSide projectionSideMatrix, ProjectionSide projectionSidePreconditioner,
            LagrangeMultiplierSeparation lagrangeSeparation)
        {
            this.maxIterationsProvider = maxIterationsProvider;
            this.pcgConvergenceTolerance = pcgConvergenceTolerance;
            this.pcgConvergenceStrategyFactory = pcgConvergenceStrategyFactory;
            this.projectionSideMatrix = projectionSideMatrix;
            this.projectionSidePreconditioner = projectionSidePreconditioner;
            this.lagrangeSeparation = lagrangeSeparation;
        }

        public Vector CalcLagrangeMultipliers(Feti1FlexibilityMatrix flexibility, IFetiPreconditioner preconditioner, 
            Feti1Projection projection, Vector disconnectedDisplacements, Vector rigidBodyModesWork, double globalForcesNorm,
            DualSolverLogger logger)
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
            pcgBuilder.MaxIterationsProvider = maxIterationsProvider;
            pcgBuilder.ResidualTolerance = pcgConvergenceTolerance;
            pcgBuilder.Convergence = pcgConvergenceStrategyFactory.CreateConvergenceStrategy(globalForcesNorm);
            PcgAlgorithm pcg = pcgBuilder.Build();
            var lagrangesBar = Vector.CreateZero(systemOrder);
            IterativeStatistics stats = pcg.Solve(pcgMatrix, pcgPreconditioner, pcgRhs, lagrangesBar, true,
                () => Vector.CreateZero(systemOrder));

            // Log statistics about PCG execution
            if (!stats.HasConverged)
            {
                throw new IterativeSolverNotConvergedException(Feti1Solver.name + " did not converge to a solution. PCG"
                    + $" algorithm run for {stats.NumIterationsRequired} iterations and the residual norm ratio was"
                    + $" {stats.ResidualNormRatioEstimation}");
            }
            logger.PcgIterations = stats.NumIterationsRequired;
            logger.PcgResidualNormRatio = stats.ResidualNormRatioEstimation;

            // Calculate the actual lagrange multipliers from the separation formula: λ = λ0 + P * λbar
            var lagranges = Vector.CreateZero(systemOrder);
            projection.ProjectVector(lagrangesBar, lagranges, false);
            lagranges.AddIntoThis(lagrangesParticular);
            return lagranges;
        }

        /// <summary>
        /// Calculates the actual lagrange multipliers depending on the separation formula.
        /// </summary>
        internal Vector CombineLagrangeMultipliers(Vector lagrangesParticular, Vector lagrangesBar, Feti1Projection projection)
        {
            if (lagrangeSeparation == LagrangeMultiplierSeparation.Simple)
            {
                // λ = λ0 + P * λbar
                return lagrangesBar + lagrangesParticular;
            }
            else if (lagrangeSeparation == LagrangeMultiplierSeparation.WithProjection)
            {
                // λ = λ0 + P * λbar
                var lagranges = Vector.CreateZero(lagrangesBar.Length);
                projection.ProjectVector(lagrangesBar, lagranges, false);
                lagranges.AddIntoThis(lagrangesParticular);
                return lagranges;
            }
            else throw new ArgumentException("The lagrange separation can only be: a) λ = λ0 + λbar or b) λ = λ0 + P * λbar");
        }

        private PcgMatrix DefinePcgMatrix(Feti1FlexibilityMatrix flexibility, Feti1Projection projection)
        {
            //TODO: projection.ProjectVector() overwrites the rhs vector passed in. However that is not clear by the method.
            //      I should either always return the resulting vector (negating any optimizations) or set the convention that
            //      vectors passed in as results will always be overwritten. The latter forces the method that accepts them to 
            //      clear them (if the method needs to), which is not needed if the resulting vector has just be initialized to 0.

            if (projectionSideMatrix == ProjectionSide.Left)
            {
                // A * x = (P^T * F) * x = P^T * (F * x)
                return new PcgMatrix(flexibility.Order, (x, y) => projection.ProjectVector(flexibility.Multiply(x), y, true));
            }
            else if (projectionSideMatrix == ProjectionSide.Both)
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
            if (projectionSidePreconditioner == ProjectionSide.None)
            {
                // x = prec(A) * y => x = prec(F) * y
                return new PcgPreconditioner((x, y) => preconditioner.SolveLinearSystem(y, x));
            }
            else if (projectionSidePreconditioner == ProjectionSide.Left)
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
            else if (projectionSidePreconditioner == ProjectionSide.Both)
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

        public class Builder
        {
            public LagrangeMultiplierSeparation LagrangeSeparation { get; set; } = LagrangeMultiplierSeparation.WithProjection;
            public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);
            public IFetiPcgConvergenceFactory PcgConvergenceStrategyFactory { get; set; } =
                new ApproximateResidualConvergence.Factory();
            public double PcgConvergenceTolerance { get; set; } = 1E-7;
            public ProjectionSide ProjectionSideMatrix { get; set; } = ProjectionSide.Both;
            public ProjectionSide ProjectionSidePreconditioner { get; set; } = ProjectionSide.Both;

            public Feti1ProjectedInterfaceProblemSolver Build() => new Feti1ProjectedInterfaceProblemSolver(
                MaxIterationsProvider, PcgConvergenceTolerance, PcgConvergenceStrategyFactory,  
                ProjectionSideMatrix, ProjectionSidePreconditioner, LagrangeSeparation);
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
