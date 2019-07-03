using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;

// TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip can be found
namespace ISAAR.MSolve.XFEM.Analyzers
{
    /// <summary>
    /// Implements crack propagation under static loading with linear material behavior. Based on Linear Elastic Fracture 
    /// Mechanics. Appropriate for brittle materials or fatigue crack propagation analysis. For now, it only works with XFEM.
    /// </summary>
    public class QuasiStaticCrackPropagationAnalyzer //: IAnalyzer
    {
        private readonly ICrackDescription crack;
        private readonly double fractureToughness;
        private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private readonly int maxIterations;
        private readonly XModel model;
        private readonly bool reanalysis = true;

        //private readonly IStaticProvider problem; //TODO: refactor and use this instead
        private readonly ElementStructuralStiffnessProvider problem = new ElementStructuralStiffnessProvider();
        private readonly DirichletEquivalentLoadsStructural loadsAssembler; 

        private readonly ISolver solver;
        private HashSet<ISubdomain> newTipEnrichedSubdomains;

        public QuasiStaticCrackPropagationAnalyzer(XModel model, ISolver solver, /*IStaticProvider problem,*/
            ICrackDescription crack, double fractureToughness, int maxIterations)
        {
            this.model = model;
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
            //this.problem = problem;
            this.crack = crack;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;

            //TODO: Refactor problem structural and remove the next
            problem = new ElementStructuralStiffnessProvider();
            loadsAssembler = new DirichletEquivalentLoadsStructural(problem); ;
        }

        public IDomainDecompositionLogger DDLogger { get; set; }
        public CrackPropagationTermination Termination { get; private set;}

        public void Initialize(bool isFirstAnalysis = true)
        {
            // The order in which the next initializations happen is very important.
            if (isFirstAnalysis) model.ConnectDataStructures();

            //solver.Initialize(); //TODO: not sure about this one.
        }

        /// <summary>
        /// Returns the crack path after repeatedly executing: XFEM analysis, SIF calculation, crack propagation
        /// </summary>
        /// <returns></returns>
        public void Analyze()
        {
            int analysisStep;
            for (analysisStep = 0; analysisStep < maxIterations; ++analysisStep)
            {
                Debug.WriteLine($"Crack propagation step {analysisStep}");
                Console.WriteLine($"Crack propagation step {analysisStep}");

                // Apply the updated enrichements.
                crack.UpdateEnrichments();

                // Update the mesh partitioning and identify unmodified subdomains to avoid fully processing them again.
                UpdateSubdomains();

                // Order and count dofs
                solver.OrderDofs(false);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    if (linearSystem.Subdomain.ConnectivityModified)
                    {
                        linearSystem.Reset(); // Necessary to define the linear system's size 
                        linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                    }
                }

                // Create the stiffness matrix and then the forces vector
                //problem.ClearMatrices();
                BuildMatrices();
                model.AssignLoads(solver.DistributeNodalLoads);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
                }
                AddEquivalentNodalLoadsToRhs();

                // Plot domain decomposition data, if necessary
                if (DDLogger != null) DDLogger.PlotSubdomains(model);

                // Solve the linear system
                solver.Solve();

                //// Output field data
                //if (fieldOutput != null)
                //{
                //    fieldOutput.WriteOutputData(solver.DofOrderer, freeDisplacements, constrainedDisplacements, iteration);
                //}

                // Let the crack propagate
                //Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(solver.DofOrderer);
                var freeDisplacements = new Dictionary<int, Vector>();
                foreach (int s in linearSystems.Keys) freeDisplacements[s] = (Vector)(linearSystems[s].Solution); //TODO: avoid this cast.
                crack.Propagate(freeDisplacements);

                // Check convergence 
                //TODO: Perhaps this should be done by the crack geometry or the Propagator itself and handled via exceptions 

                foreach (var tipPropagator in crack.CrackTipPropagators)
                {
                    double sifEffective = CalculateEquivalentSIF(tipPropagator.Value.Logger.SIFsMode1[analysisStep],
                        tipPropagator.Value.Logger.SIFsMode2[analysisStep]);
                    //Console.WriteLine("Keff = " + sifEffective);
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

        private void AddEquivalentNodalLoadsToRhs()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                try
                {
                    // Make sure there is at least one non zero prescribed displacement.
                    (INode node, IDofType dof, double displacement) = linearSystem.Subdomain.Constraints.Find(du => du != 0.0);

                    //TODO: the following 2 lines are meaningless outside diplacement control (and even then, they are not so clear).
                    double scalingFactor = 1;
                    IVector initialFreeSolution = linearSystem.CreateZeroVector();

                    //IVector equivalentNodalLoads = problem.DirichletLoadsAssembler.GetEquivalentNodalLoads(
                    //    linearSystem.Subdomain, initialFreeSolution, scalingFactor);
                    IVector equivalentNodalLoads = loadsAssembler.GetEquivalentNodalLoads(linearSystem.Subdomain, 
                        initialFreeSolution, scalingFactor);
                    linearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads);
                }
                catch (KeyNotFoundException)
                {
                    // There aren't any non zero prescribed displacements, therefore we do not have to calculate the equivalent 
                    // nodal loads, which is an expensive operation (all elements are accessed, their stiffness is built, etc..)
                }
            }
        }

        private void BuildMatrices()
        {
            //Dictionary<int, IMatrixView> stiffnesses = problem.CalculateMatrix();
            //foreach (ILinearSystem linearSystem in linearSystems.Values)
            //{
            //    linearSystem.Matrix = stiffnesses[linearSystem.Subdomain.ID];
            //}

            Dictionary<int, IMatrix> stiffnesses = solver.BuildGlobalMatrices(problem);
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.Matrix = stiffnesses[linearSystem.Subdomain.ID];
            }
        }

        // TODO: Abstract this and add Tanaka_1974 approach
        private double CalculateEquivalentSIF(double sifMode1, double sifMode2)
        {
            return Math.Sqrt(sifMode1 * sifMode1 + sifMode2 * sifMode2);
        }

        private HashSet<ISubdomain> FindSubdomainsWithNewHeavisideEnrichedNodes()
        {
            var newHeavisideEnrichedSubdomains = new HashSet<ISubdomain>();
            foreach (ISet<XNode> heavisideNodes in crack.CrackBodyNodesNew.Values)
            {
                foreach (XNode node in heavisideNodes)
                {
                    newHeavisideEnrichedSubdomains.UnionWith(node.SubdomainsDictionary.Values);
                }
            }
            return newHeavisideEnrichedSubdomains;
        }

        private HashSet<ISubdomain> FindSubdomainsWithNewTipEnrichedNodes()
        {
            var newTipEnrichedSubdomains = new HashSet<ISubdomain>();
            foreach (ISet<XNode> tipNodes in crack.CrackTipNodesNew.Values)
            {
                foreach (XNode node in tipNodes)
                {
                    newTipEnrichedSubdomains.UnionWith(node.SubdomainsDictionary.Values);
                }
            }
            return newTipEnrichedSubdomains;
        }

        private void UpdateSubdomains()
        {
            if (model.Subdomains.Count == 1) return;

            if (newTipEnrichedSubdomains == null) 
            {
                // First analysis step: All subdomains must be fully processed.
                foreach (XSubdomain subdomain in model.Subdomains.Values)
                {
                    subdomain.ConnectivityModified = true;
                    subdomain.StiffnessModified = true;
                }

                // Prepare for the next analysis step
                newTipEnrichedSubdomains = FindSubdomainsWithNewTipEnrichedNodes();
            }
            else
            {
                // Update the mesh partitioning, if necessary
                HashSet<ISubdomain> modifiedSubdomains;
                modifiedSubdomains = new HashSet<ISubdomain>();

                if (reanalysis)
                {
                    // Set them all to unmodified.
                    foreach (XSubdomain subdomain in model.Subdomains.Values)
                    {
                        subdomain.ConnectivityModified = false;
                        subdomain.StiffnessModified = false;
                    }

                    // The modified subdomains are the ones containing nodes enriched with tip or Heaviside functions during the  
                    // current analysis step. Also the ones that had tip enriched nodes in the previous step.
                    modifiedSubdomains.UnionWith(newTipEnrichedSubdomains);
                    newTipEnrichedSubdomains = FindSubdomainsWithNewTipEnrichedNodes(); // Prepare for the next analysis step
                    modifiedSubdomains.UnionWith(newTipEnrichedSubdomains);
                    HashSet<ISubdomain> newHeavisideSubdomains = FindSubdomainsWithNewHeavisideEnrichedNodes();
                    modifiedSubdomains.UnionWith(newHeavisideSubdomains);

                    foreach (ISubdomain subdomain in modifiedSubdomains)
                    {
                        subdomain.ConnectivityModified = true;
                        subdomain.StiffnessModified = true;
                    }
                }
            }
        }
    }
}
