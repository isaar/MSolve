using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019
{
    public static class Utilities
    {
        public enum SolverType { Skyline, Feti1, FetiDP }

        public enum BoundaryConditions
        {
            BottomConstrainXY_TopDisplacementY, BottomConstrainXY_TopConstrainXDisplacementY,
            BottomConstrainY_TopDisplacementY, BottomDisplacementY_TopDisplacementY,
            BottomConstrainXDisplacementY_TopConstrainXDisplacementY
        }

        public static Dictionary<int, INode[]> FindCornerNodesFromCrosspoints2D(IStructuralModel model)
        {
            //TODO: This is also done by the analyzer. Perhaps it should not.
            model.ConnectDataStructures();
            var cornerNodes = new Dictionary<int, INode[]>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                // In 2D, if multiplicity > 2, the node is a crosspoint 
                cornerNodes[subdomain.ID] = subdomain.Nodes.Where(node => node.SubdomainsDictionary.Count > 2).ToArray();
            }
            return cornerNodes;
        }

        public static void PlotSubdomains(string plotPath, XModel model)
        {
            model.ConnectDataStructures();
            var writer = new MeshPartitionWriter();
            var nodesPerSubdomain = new Dictionary<int, IReadOnlyList<XNode>>();
            var elementsPerSubdomain = new Dictionary<int, IReadOnlyList<IXFiniteElement>>();
            foreach (int subdomainID in model.Subdomains.Keys)
            {
                nodesPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Nodes;
                elementsPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Elements;
            }
            writer.WriteSubdomainElements(plotPath, nodesPerSubdomain, elementsPerSubdomain);
        }

        public static double RunSingleCrackedStep(XModel model, ICrackDescription crack, ISolver solver)
        {
            // Enrich nodes to take the crack into account
            crack.UpdateEnrichments();

            var problem = new ProblemStructural(model, solver);
            var linearAnalyzer = new LinearAnalyzer(model, solver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, solver, problem, linearAnalyzer);

            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Return norm2(displacements)
            if (model.Subdomains.Count == 1)
            {
                return solver.LinearSystems.First().Value.Solution.Norm2();
            }
            else if (solver is Feti1Solver feti1Solver)
            {
                var sudomainDisplacements = new Dictionary<int, IVectorView>();
                foreach (var ls in feti1Solver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
                return feti1Solver.GatherGlobalDisplacements(sudomainDisplacements).Norm2();
            }
            else if (solver is FetiDPSolver fetiDPSolver)
            {
                var sudomainDisplacements = new Dictionary<int, IVectorView>();
                foreach (var ls in fetiDPSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
                return fetiDPSolver.GatherGlobalDisplacements(sudomainDisplacements).Norm2();
            }
            else throw new NotImplementedException("Invalid solver");
        }

        public static double RunUncrackedAnalysis(XModel model, ISolver solver)
        {
            var problem = new ProblemStructural(model, solver);
            var linearAnalyzer = new LinearAnalyzer(model, solver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, solver, problem, linearAnalyzer);

            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Return norm2(displacements)
            if (model.Subdomains.Count == 1)
            {
                return solver.LinearSystems.First().Value.Solution.Norm2();
            }
            else if (solver is Feti1Solver feti1Solver)
            {
                var sudomainDisplacements = new Dictionary<int, IVectorView>();
                foreach (var ls in feti1Solver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
                return feti1Solver.GatherGlobalDisplacements(sudomainDisplacements).Norm2();
            }
            else if (solver is FetiDPSolver fetiDPSolver)
            {
                var sudomainDisplacements = new Dictionary<int, IVectorView>();
                foreach (var ls in fetiDPSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
                return fetiDPSolver.GatherGlobalDisplacements(sudomainDisplacements).Norm2();
            }
            else throw new NotImplementedException("Invalid solver");
        }
    }
}
