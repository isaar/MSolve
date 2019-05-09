using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// base class for Primary multiscale analysis classes that connect nesessary structures for a FE2 simulation
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public abstract class StructuralProblemsMicrostructureBase
    {
        public int SolverData { get; set; }

        public virtual Dictionary<int,Element> GetBoundaryFiniteElementsDictionary(Model model, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Element> boundaryElements = new Dictionary<int, Element>();

            foreach(Element element in model.Elements)
            {
                bool containsBoundaryNode = false;

                var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
                var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
                for (int j = 0; j < elementDOFTypes.Count; j++)
                {
                    INode elementNode = matrixAssemblyNodes[j];
                    if (boundaryNodes.ContainsKey(elementNode.ID))
                    {
                        containsBoundaryNode = true;
                        break;
                    }
                }

                if (containsBoundaryNode)
                {
                    boundaryElements.Add(element.ID, element);
                }

            }

            return boundaryElements;

        }

        public virtual Dictionary<int, Element> GetBoundaryFiniteElementsDictionary(Subdomain subdomain, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Element> subdomainBoundaryElements = new Dictionary<int, Element>();

            foreach (Element element in subdomain.Elements)
            {
                bool containsBoundaryNode = false;

                var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
                var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
                for (int j = 0; j < elementDOFTypes.Count; j++)
                {
                    INode elementNode = matrixAssemblyNodes[j];
                    if (boundaryNodes.ContainsKey(elementNode.ID))
                    {
                        containsBoundaryNode = true;
                        break;
                    }
                }

                if (containsBoundaryNode)
                {
                    subdomainBoundaryElements.Add(element.ID, element);
                }

            }

            return subdomainBoundaryElements;

        }

        public virtual Dictionary<int, Dictionary<int, Element>> GetSubdomainsBoundaryFiniteElementsDictionaries(Model model, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Dictionary<int, Element>> subdomainsBoundaryElements = new Dictionary<int, Dictionary<int, Element>>();

            foreach (Subdomain subdomain in model.Subdomains)
            {
                Dictionary<int, Element> subdBoundaryElements = GetBoundaryFiniteElementsDictionary(subdomain, boundaryNodes);
                subdomainsBoundaryElements.Add(subdomain.ID, subdBoundaryElements);
            }

            return subdomainsBoundaryElements;
            
        }

        public virtual (MicrostructureBvpNRNLAnalyzer, ProblemStructural,ElementStructuralStiffnessProvider) AnalyzeMicrostructure(Model model,  ISolver solver,
            int increments, int MaxIterations, int IterationsForMatrixRebuild, Dictionary<int, Dictionary<IDofType, double>> totalPrescribedBoundaryDisplacements,
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Node> boundaryNodes, Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain)
        {
            IReadOnlyDictionary<int, ILinearSystem> linearSystems = solver.LinearSystems; //V2.1

            #region Creation of nessesary analyzers for NRNLAnalyzer
            ProblemStructural provider = new ProblemStructural(model, solver);

            var subdomainUpdaters = new Dictionary<int, NonLinearSubdomainUpdaterWithInitialConditions>(1); //v2.2
            //var subdomainUpdaters = new NonLinearSubdomainUpdaterWithInitialConditions[totalSubdomains];
            
            foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)"
            {
                subdomainUpdaters.Add(subdomain.ID, new NonLinearSubdomainUpdaterWithInitialConditions(subdomain)); //v2.3
                //subdomainUpdaters[counter] = new NonLinearSubdomainUpdaterWithInitialConditions(subdomain);
            }

            ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();

            //v2.4
            Dictionary<int, EquivalentContributionsAssebler> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            //equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                equivalentContributionsAssemblers.Add(subdomain.ID, new EquivalentContributionsAssebler(subdomain, elementProvider)); //v2.5
            }
            #endregion

            #region Creation of Microstructure analyzer (NRNLdevelop temporarilly). 
            MicrostructureBvpNRNLAnalyzer microAnalyzer = new MicrostructureBvpNRNLAnalyzer(model,solver, subdomainUpdaters, 
                provider, increments,  uInitialFreeDOFDisplacementsPerSubdomain,
                boundaryNodes, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, equivalentContributionsAssemblers);
            microAnalyzer.SetMaxIterations = MaxIterations;
            microAnalyzer.SetIterationsForMatrixRebuild = IterationsForMatrixRebuild;
            #endregion

            #region solution and update ------------->THA MPEI ENTOS KLASHS: of free converged displacements vectors;
            MSParentAnalyzer parentAnalyzer = new MSParentAnalyzer(model, solver, provider, microAnalyzer);
            //parentAnalyzer.BuildMatrices(); //v2.6 ston neon static analyzer den xreiazetai to build matrices poia
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion

            return (microAnalyzer,provider,elementProvider);
        }
    }
}
