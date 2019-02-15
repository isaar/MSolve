using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// base class for Primary multiscale analysis classes that connect nesessary structures for a FE2 simulation
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public abstract class StructuralProblemsMicrostructureBase_v2
    {
        public int SolverData { get; set; }

        #region v1 will not be needed in v2
        public virtual Dictionary<int, ILinearSystem> CreateNecessaryLinearSystems(Model model)
        {
            var linearSystems = new Dictionary<int, ILinearSystem>();
            foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)" tou opoiu ta stoixeia ginontai access kai me ID
            {
                linearSystems.Add(subdomain.ID, new SkylineLinearSystem(subdomain.ID, subdomain.Forces));// prosoxh sto Id twn subdomain
            }

            return linearSystems;
        }

        public virtual ISolver GetAppropriateSolver(Dictionary<int, ILinearSystem> linearSystems)
        {
            //TODO: use solver data to create the chosen ISolver
            if (linearSystems.Keys.Count == 1)
            {
                var solver = new SolverSkyline(linearSystems[1]); //linearSystems.ElementAt(0);
                return solver;
            }
            else
            {
                throw new NotImplementedException();
            }
        }
        #endregion

        public virtual Dictionary<int,Element_v2> GetBoundaryFiniteElementsDictionary_v2(Model_v2 model, Dictionary<int, Node_v2> boundaryNodes)
        {
            Dictionary<int, Element_v2> boundaryElements = new Dictionary<int, Element_v2>();

            foreach(Element_v2 element in model.Elements)
            {
                bool containsBoundaryNode = false;

                var elementDOFTypes = element.ElementType.DofEnumerator.GetDOFTypes(element);
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

        public virtual Dictionary<int, Element_v2> GetBoundaryFiniteElementsDictionary(Subdomain_v2 subdomain, Dictionary<int, Node_v2> boundaryNodes)
        {
            Dictionary<int, Element_v2> subdomainBoundaryElements = new Dictionary<int, Element_v2>();

            foreach (Element_v2 element in subdomain.Elements)
            {
                bool containsBoundaryNode = false;

                var elementDOFTypes = element.ElementType.DofEnumerator.GetDOFTypes(element);
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

        public virtual Dictionary<int, Dictionary<int, Element_v2>> GetSubdomainsBoundaryFiniteElementsDictionaries_v2(Model_v2 model, Dictionary<int, Node_v2> boundaryNodes)
        {
            Dictionary<int, Dictionary<int, Element_v2>> subdomainsBoundaryElements = new Dictionary<int, Dictionary<int, Element_v2>>();

            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                Dictionary<int, Element_v2> subdBoundaryElements = GetBoundaryFiniteElementsDictionary(subdomain, boundaryNodes);
                subdomainsBoundaryElements.Add(subdomain.ID, subdBoundaryElements);
            }

            return subdomainsBoundaryElements;
            
        }

        public virtual (MicrostructureBvpNRNLAnalyzer, ProblemStructural_v2,ElementStructuralStiffnessProvider_v2) AnalyzeMicrostructure_v2(Model_v2 model,  ISolver_v2 solver,
            int increments, int MaxIterations, int IterationsForMatrixRebuild, Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Node_v2> boundaryNodes, Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain)
        {
            IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems = solver.LinearSystems; //V2.1

            #region Creation of nessesary analyzers for NRNLAnalyzer
            ProblemStructural_v2 provider = new ProblemStructural_v2(model, solver);

            var subdomainUpdaters = new Dictionary<int, NonLinearSubdomainUpdaterWithInitialConditions_v2>(1); //v2.2
            //var subdomainUpdaters = new NonLinearSubdomainUpdaterWithInitialConditions[totalSubdomains];
            
            foreach (Subdomain_v2 subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)"
            {
                subdomainUpdaters.Add(subdomain.ID, new NonLinearSubdomainUpdaterWithInitialConditions_v2(subdomain)); //v2.3
                //subdomainUpdaters[counter] = new NonLinearSubdomainUpdaterWithInitialConditions(subdomain);
            }

            ElementStructuralStiffnessProvider_v2 elementProvider = new ElementStructuralStiffnessProvider_v2();

            //v2.4
            Dictionary<int, EquivalentContributionsAssebler_v2> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler_v2>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            //equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));
            foreach (Subdomain_v2 subdomain in model.SubdomainsDictionary.Values)
            {
                equivalentContributionsAssemblers.Add(subdomain.ID, new EquivalentContributionsAssebler_v2(subdomain, elementProvider)); //v2.5
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
