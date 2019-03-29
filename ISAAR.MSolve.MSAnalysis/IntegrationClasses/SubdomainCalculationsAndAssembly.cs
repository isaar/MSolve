using System;
using System.Collections.Generic;
//using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.Interfaces;
using System.Linq;
using ISAAR.MSolve.FEM;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysisMerge;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Supportive class  that implements nesessary integration methods associated with FE2 multiscale analysis 
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class SubdomainCalculationsAndAssembly
    {
        //SubdomainCalculationsSimultaneousObje_v2

        private Dictionary<int, double[][]> KfpDqVectors;
        private Dictionary<int, double[][]> KppDqVectors;
        //ISubdomainFreeDofOrdering dofOrdering;
        //DofTable FreeDofs;
        //v2.1 Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary;
        IScaleTransitions_v2 scaleTransitions;
        Dictionary<int, Dictionary<int, Element_v2>> boundaryElements;
        Dictionary<int, Node_v2> boundaryNodes;
        Dictionary<int, int> boundaryNodesOrder;
        //int currentSubdomainID;

        public (Dictionary<int, double[][]>, Dictionary<int, double[][]>) UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultipleObje_v2(Model_v2 model, IElementMatrixProvider_v2 elementProvider, IScaleTransitions_v2 scaleTransitions,
            Dictionary<int, Node_v2> boundaryNodes, Dictionary<int, Dictionary<int, Element_v2>> boundaryElements,ISolver_v2 solver)
        {
            IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems = solver.LinearSystems; //v2.3

            Dictionary<int, double[][]> KfpDqSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            Dictionary<int, double[][]> KppDqVectorsSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            this.boundaryElements = boundaryElements;
            this.boundaryNodes = boundaryNodes;
            this.scaleTransitions = scaleTransitions;

            KfpDqVectors = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            KppDqVectors = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                #region Create KfpDq and KppDq vectors 
                KfpDqVectors[subdomain.ID] = new double[scaleTransitions.MacroscaleVariableDimension()][];
                for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
                {
                    KfpDqVectors[subdomain.ID][j1] = new double[subdomain.FreeDofOrdering.NumFreeDofs]; //v2.2 subdomain.TotalDOFs]; 
                }

                KppDqVectors[subdomain.ID] = new double[scaleTransitions.MacroscaleVariableDimension()][];
                boundaryNodesOrder = SubdomainCalculations_v2.GetNodesOrderInDictionary(boundaryNodes);
                for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
                {
                    KppDqVectors[subdomain.ID][j1] = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)
                }
                #endregion
            }

            var StiffnessProvider = new StiffnessProviderSimu_v2(this);
            Dictionary<int, IMatrix> subdomainKs = solver.BuildGlobalMatrices(StiffnessProvider);

            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                //dofOrdering = subdomain.FreeDofOrdering; //_v2.1
                //FreeDofs = subdomain.FreeDofOrdering.FreeDofs;//_v2.1 nodalDOFsDictionary = subdomain.NodalDOFsDictionary;
                //currentSubdomainID = subdomain.ID;

                


                //v2.4 var subdomainK= GlobalMatrixAssemblerSkyline.CalculateFreeFreeGlobalMatrix(subdomain, StiffnessProvider);                

                linearSystems[subdomain.ID].Matrix = subdomainKs[subdomain.ID];
                //v2.5 linearSystems[subdomain.ID].Matrix = subdomainK;

                KfpDqSubdomains.Add(subdomain.ID, KfpDqVectors[subdomain.ID]);
                KppDqVectorsSubdomains.Add(subdomain.ID, KppDqVectors[subdomain.ID]);                
            }

            return (KfpDqSubdomains, KppDqVectorsSubdomains);        
        }

        public void UpdateVectors_v2(IElement_v2 element, IMatrix ElementK)
        {
            ISubdomain_v2 subdomain = element.Subdomain;
            if (boundaryElements[subdomain.ID].ContainsKey(element.ID))//COPIED From UpdateSubdomainKffAndCalculateKfpDqAndKppDqp (prosoxh boundary elements Dictionary diathetoun kai to model kai to subdomain kai einai diaforetika edw exei diorthwthei
            {
                //ADDED these lines from another part of UpdateSubdomainKffAndCalculateKfpDqAndKppDqp
                var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
                var elementDOFTypes = element.ElementType.DofEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);


                #region KfpDq Multiplication
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    int dofTypeRowToNumber = -1; //v2.6
                    foreach (DOFType dofTypeRow in elementDOFTypes[i])
                    {
                        dofTypeRowToNumber++;
                        bool isFree = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(matrixAssemblyNodes[i], elementDOFTypes[i][dofTypeRowToNumber],
                        out int dofRow); //v2.6
                        //int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (isFree) // TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                        {                    // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                //foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                //{
                                //    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                //    if (dofColumn != -1)
                                //    {
                                //        int height = dofRow - dofColumn;
                                //        if (height >= 0)
                                //            K.Data[K.RowIndex[dofRow] + height] += ElementK[iElementMatrixRow, iElementMatrixColumn];
                                //    }
                                //    iElementMatrixColumn++;
                                //}
                                int nodalDofsNumber = elementDOFTypes[j].Count; //TODOGerasimos elegxos oti edw oi ginetai prosvash apo 0:1:megethos 
                                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                                {
                                    double[] element_Kfp_triplette = new double[nodalDofsNumber]; //nodalDofsNumber: giati oxi scaleTransitions.PrescribedDofsPerNode()? Dioti tou ta pairname ola(triplette) kai dialegei to 
                                    for (int j1 = 0; j1 < nodalDofsNumber; j1++)                  //scaleTransitions.MicroToMacroTransition ti tha xrhsimopoihsei apo afta analoga pws einai implemented
                                    {
                                        element_Kfp_triplette[j1] = ElementK[iElementMatrixRow, iElementMatrixColumn + j1];
                                    }

                                    double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kfp_triplette);
                                    for (int j2 = 0; j2 < contribution.GetLength(0); j2++)
                                    {
                                        KfpDqVectors[subdomain.ID][j2][dofRow] += contribution[j2]; // TODO diorthothike
                                    }

                                }
                                iElementMatrixColumn += nodalDofsNumber;

                            }
                        }
                        iElementMatrixRow++;
                    }
                }
                #endregion

                #region KppDq Multiplications
                iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodes.ContainsKey(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                //
                                //foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                //{
                                //    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                //    if (dofColumn != -1)// TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                                //    {                   // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok

                                //        for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
                                //        {
                                //            f3_vectors[i2][dofrow_p] += ElementK[iElementMatrixRow + i1, iElementMatrixColumn] * f2_vectors[i2][dofColumn]; /////
                                //        }

                                //    }
                                //    iElementMatrixColumn++;
                                //}
                                //
                                int nodalDofsNumber = elementDOFTypes[j].Count;
                                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                                {
                                    double[] element_Kpp_triplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                                    for (int j2 = 0; j2 < scaleTransitions.PrescribedDofsPerNode(); j2++)
                                    {
                                        element_Kpp_triplette[j2] = ElementK[iElementMatrixRow + i1, iElementMatrixColumn + j2]; // mallon iElementMatrixRow + i1
                                    }
                                    double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kpp_triplette);
                                    for (int j1 = 0; j1 < contribution.GetLength(0); j1++)
                                    {
                                        KppDqVectors[subdomain.ID][j1][dofrow_p] += contribution[j1];
                                    }
                                }
                                iElementMatrixColumn += nodalDofsNumber;



                            }

                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                #endregion
            }
        }

    }
}
