using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;

// TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip can be found
namespace ISAAR.MSolve.Analyzers.CrackPropagation
{
    /// <summary>
    /// Implements crack propagation under static loading with linear material behavior. Based on Linear Elastic Fracture 
    /// Mechanics. Appropriate for brittle materials or fatigue crack propagation analysis. For now, it only works with XFEM.
    /// </summary>
    public class QuasiStaticCrackPropagationAnalyzer : IAnalyzer
    {
        private readonly ICrackDescription crack;
        private readonly double fractureToughness;
        private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private readonly int maxIterations;
        private readonly XModel model;
        private readonly IStaticProvider provider;
        private readonly ISolver solver;

        public QuasiStaticCrackPropagationAnalyzer(XModel model, ISolver solver, IStaticProvider provider,
            ICrackDescription crack, double fractureToughness, int maxIterations)
        {
            this.model = model;
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
            this.provider = provider;
            this.crack = crack;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
        }

        public Dictionary<int, IAnalyzerLog[]> Logs => throw new NotImplementedException();

        public CrackPropagationTermination Termination { get; private set;}

        public void BuildMatrices()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.Matrix = provider.CalculateMatrix(linearSystem.Subdomain);
            }
        }

        public void Initialize(bool isFirstAnalysis = true)
        {
            // The order in which the next initializations happen is very important.
            if (isFirstAnalysis) model.ConnectDataStructures();

            solver.Initialize(); //TODO: not sure about this one.
        }

        /// <summary>
        /// Returns the crack path after repeatedly executing: XFEM analysis, SIF calculation, crack propagation
        /// </summary>
        /// <returns></returns>
        public void Analyze()
        {
            int iteration;
            for (iteration = 0; iteration < maxIterations; ++iteration)
            {
                // Apply the updated enrichements.
                crack.UpdateEnrichments();

                // Order and count dofs
                solver.OrderDofs(false);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }

                // Create the stiffness matrix and then the forces vector
                BuildMatrices();
                model.AssignLoads(solver.DistributeNodalLoads);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
                }

                // Solve the linear system
                solver.Solve();
                //Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(solver.DofOrderer);
                Vector freeDisplacements = (Vector)(linearSystems[0].Solution); //TODO: avoid this.

                //// Output field data
                //if (fieldOutput != null)
                //{
                //    fieldOutput.WriteOutputData(solver.DofOrderer, freeDisplacements, constrainedDisplacements, iteration);
                //}

                // Let the crack propagate
                crack.Propagate(freeDisplacements);

                // Check convergence 
                //TODO: Perhaps this should be done by the crack geometry or the Propagator itself and handled via exceptions 

                foreach (var tipPropagator in crack.CrackTipPropagators)
                {
                    double sifEffective = CalculateEquivalentSIF(tipPropagator.Value.Logger.SIFsMode1[iteration],
                    tipPropagator.Value.Logger.SIFsMode2[iteration]);
                    if (sifEffective >= fractureToughness)
                    {
                        Termination = CrackPropagationTermination.FractureToughnessIsExceeded;
                        return;
                    }
                    if (!model.Boundary.IsInside(tipPropagator.Key))
                    {
                        Termination = CrackPropagationTermination.CrackExitsDomainBoundary;
                        return;
                    }
                }
            }
            Termination = CrackPropagationTermination.RequiredIterationsWereCompleted;
        }

        // TODO: Abstract this and add Tanaka_1974 approach
        private double CalculateEquivalentSIF(double sifMode1, double sifMode2)
        {
            return Math.Sqrt(sifMode1 * sifMode1 + sifMode2 * sifMode2);
        }

        public void Solve()
        {
            // TODO: This is necessary for XFEM, but other approaches may not need to reorder dofs at each iteration.
            solver.OrderDofs(false);
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.Reset(); // Necessary to define the linear system's size 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
            }

            BuildMatrices();

            // Loads must be created after building the matrices.
            //TODO: Some loads may not have to be recalculated each time the stiffness changes.
            model.AssignLoads(solver.DistributeNodalLoads);
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }
        }
    }
}
