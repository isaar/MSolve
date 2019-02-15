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
    public class SubdomainCalculationsSimultaneousObje_v2
    {
        //public static (double[][], double[][], SkylineMatrix2D) UpdateSubdomainKffAndCalculateKfpDqAndKppDqp(Subdomain subdomain, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions,
        //    Dictionary<int, Node> boundaryNodes, Dictionary<int, Element> boundaryElements)
        //{
        //    var nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

        //    #region Create KfpDq and KppDq vectors 
        //    double[][] KfpDqVectors = new double[scaleTransitions.MacroscaleVariableDimension()][];
        //    for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
        //    {
        //        KfpDqVectors[j1] = new double[subdomain.TotalDOFs]; // h allliws subdomain.Forces.GetLength(0)
        //    }

        //    double[][] KppDqVectors = new double[scaleTransitions.MacroscaleVariableDimension()][];
        //    Dictionary<int, int> boundaryNodesOrder = SubdomainCalculations_v2.GetNodesOrderInDictionary(boundaryNodes);
        //    for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
        //    {
        //        KppDqVectors[j1] = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)
        //    }
        //    #endregion


        //    var times = new Dictionary<string, TimeSpan>();
        //    var totalStart = DateTime.Now;
        //    SkylineMatrix2D K = new SkylineMatrix2D(CalculateRowIndex(subdomain, nodalDOFsDictionary));
        //    times.Add("rowIndexCalculation", DateTime.Now - totalStart);
        //    times.Add("element", TimeSpan.Zero);
        //    times.Add("addition", TimeSpan.Zero);
        //    foreach (Element element in subdomain.ElementsDictionary.Values)
        //    {
        //        var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
        //        var elStart = DateTime.Now;
        //        IMatrix2D ElementK = elementProvider.Matrix(element);
        //        times["element"] += DateTime.Now - elStart;

        //        elStart = DateTime.Now;
        //        var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
        //        var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);

        //        #region kFF Assembly
        //        int iElementMatrixRow = 0;
        //        for (int i = 0; i < elementDOFTypes.Count; i++)
        //        {
        //            INode nodeRow = matrixAssemblyNodes[i];
        //            foreach (DOFType dofTypeRow in elementDOFTypes[i])
        //            {
        //                int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
        //                if (dofRow != -1)
        //                {
        //                    int iElementMatrixColumn = 0;
        //                    for (int j = 0; j < elementDOFTypes.Count; j++)
        //                    {
        //                        INode nodeColumn = matrixAssemblyNodes[j];
        //                        foreach (DOFType dofTypeColumn in elementDOFTypes[j])
        //                        {
        //                            int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
        //                            if (dofColumn != -1)
        //                            {
        //                                int height = dofRow - dofColumn;
        //                                if (height >= 0)
        //                                    K.Data[K.RowIndex[dofRow] + height] += ElementK[iElementMatrixRow, iElementMatrixColumn];
        //                            }
        //                            iElementMatrixColumn++;
        //                        }
        //                    }
        //                }
        //                iElementMatrixRow++;
        //            }
        //        }
        //        #endregion
                
        //        if(boundaryElements.ContainsKey(element.ID))
        //        {
        //            #region KfpDq Multiplication
        //            iElementMatrixRow = 0;
        //            for (int i = 0; i < elementDOFTypes.Count; i++)
        //            {
        //                INode nodeRow = matrixAssemblyNodes[i];
        //                foreach (DOFType dofTypeRow in elementDOFTypes[i])
        //                {
        //                    int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
        //                    if (dofRow != -1) // TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
        //                    {                    // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
        //                        int iElementMatrixColumn = 0;
        //                        for (int j = 0; j < elementDOFTypes.Count; j++)
        //                        {
        //                            INode nodeColumn = matrixAssemblyNodes[j];
        //                            //foreach (DOFType dofTypeColumn in elementDOFTypes[j])
        //                            //{
        //                            //    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
        //                            //    if (dofColumn != -1)
        //                            //    {
        //                            //        int height = dofRow - dofColumn;
        //                            //        if (height >= 0)
        //                            //            K.Data[K.RowIndex[dofRow] + height] += ElementK[iElementMatrixRow, iElementMatrixColumn];
        //                            //    }
        //                            //    iElementMatrixColumn++;
        //                            //}
        //                            int nodalDofsNumber = elementDOFTypes[j].Count; //TODOGerasimos elegxos oti edw oi ginetai prosvash apo 0:1:megethos 
        //                            if (boundaryNodes.ContainsKey(nodeColumn.ID))
        //                            {
        //                                double[] element_Kfp_triplette = new double[nodalDofsNumber]; //nodalDofsNumber: giati oxi scaleTransitions.PrescribedDofsPerNode()? Dioti tou ta pairname ola(triplette) kai dialegei to 
        //                                for (int j1 = 0; j1 < nodalDofsNumber; j1++)                  //scaleTransitions.MicroToMacroTransition ti tha xrhsimopoihsei apo afta analoga pws einai implemented
        //                                {
        //                                    element_Kfp_triplette[j1] = ElementK[iElementMatrixRow, iElementMatrixColumn + j1];
        //                                }

        //                                double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kfp_triplette);
        //                                for (int j2 = 0; j2 < contribution.GetLength(0); j2++)
        //                                {
        //                                    KfpDqVectors[j2][dofRow] += contribution[j2]; // TODO diorthothike
        //                                }

        //                            }
        //                            iElementMatrixColumn += nodalDofsNumber;

        //                        }
        //                    }
        //                    iElementMatrixRow++;
        //                }
        //            }
        //            #endregion

        //            #region KppDq Multiplications
        //            iElementMatrixRow = 0;
        //            for (int i = 0; i < elementDOFTypes.Count; i++)
        //            {
        //                INode nodeRow = matrixAssemblyNodes[i];
        //                if (boundaryNodes.ContainsKey(nodeRow.ID))
        //                {
        //                    for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
        //                    {
        //                        int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
        //                        int iElementMatrixColumn = 0;
        //                        for (int j = 0; j < elementDOFTypes.Count; j++)
        //                        {
        //                            INode nodeColumn = matrixAssemblyNodes[j];
        //                            //
        //                            //foreach (DOFType dofTypeColumn in elementDOFTypes[j])
        //                            //{
        //                            //    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
        //                            //    if (dofColumn != -1)// TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
        //                            //    {                   // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok

        //                            //        for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
        //                            //        {
        //                            //            f3_vectors[i2][dofrow_p] += ElementK[iElementMatrixRow + i1, iElementMatrixColumn] * f2_vectors[i2][dofColumn]; /////
        //                            //        }

        //                            //    }
        //                            //    iElementMatrixColumn++;
        //                            //}
        //                            //
        //                            int nodalDofsNumber = elementDOFTypes[j].Count;
        //                            if (boundaryNodes.ContainsKey(nodeColumn.ID))
        //                            {
        //                                double[] element_Kpp_triplette = new double[scaleTransitions.PrescribedDofsPerNode()];
        //                                for (int j2 = 0; j2 < scaleTransitions.PrescribedDofsPerNode(); j2++)
        //                                {
        //                                    element_Kpp_triplette[j2] = ElementK[iElementMatrixRow + i1, iElementMatrixColumn + j2]; // mallon iElementMatrixRow + i1
        //                                }
        //                                double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kpp_triplette);
        //                                for (int j1 = 0; j1 < contribution.GetLength(0); j1++)
        //                                {
        //                                    KppDqVectors[j1][dofrow_p] += contribution[j1];
        //                                }
        //                            }
        //                            iElementMatrixColumn += nodalDofsNumber;



        //                        }

        //                    }
        //                }
        //                iElementMatrixRow += elementDOFTypes[i].Count;
        //            }
        //            #endregion
        //        }

        //        times["addition"] += DateTime.Now - elStart;
        //    }
        //    var totalTime = DateTime.Now - totalStart;


        //    //linearSystems[subdomain.ID].Matrix = K;
        //    return (KfpDqVectors, KppDqVectors,K);
        //}

        //#region GlobalMatrixAssemblerSkkyline methods
        //private static int[] CalculateRowIndex(Subdomain subdomain, Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary)
        //{
        //    int[] rowHeights = new int[subdomain.TotalDOFs];
        //    foreach (Element element in subdomain.ElementsDictionary.Values)
        //    {
        //        int minDOF = Int32.MaxValue;
        //        //foreach (Node node in element.NodesDictionary.Values.Where(e => e.EmbeddedInElement == null))
        //        foreach (Node node in element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
        //        {
        //            if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement))
        //                continue;
        //            foreach (int dof in nodalDOFsDictionary[node.ID].Values)
        //                if (dof != -1) minDOF = Math.Min(dof, minDOF);
        //        }
        //        //foreach (Node node in element.NodesDictionary.Values.Where(e => e.EmbeddedInElement == null))
        //        foreach (Node node in element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
        //        {
        //            if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement))
        //                continue;
        //            foreach (int dof in nodalDOFsDictionary[node.ID].Values)
        //                if (dof != -1) rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
        //        }
        //    }

        //    int[] rowIndex = new int[subdomain.TotalDOFs + 1];
        //    rowIndex[0] = 0;
        //    rowIndex[1] = 1;
        //    for (int i = 1; i < subdomain.TotalDOFs; i++)
        //        rowIndex[i + 1] = rowIndex[i] + rowHeights[i] + 1;
        //    return rowIndex;
        //}
        //#endregion

        //public static (Dictionary<int, double[][]> , Dictionary<int, double[][]> ) UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultiple(Model model, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions,
        //    Dictionary<int, Node> boundaryNodes, Dictionary<int, Dictionary<int, Element>> boundaryElements, Dictionary<int, ILinearSystem> linearSystems)
        //{
        //    //TODO: afth h methodos antikathista (mallon) tis 227-228, 232 kai 271 grammes ths Microstru3developMultipleSubdomainsUseBase

        //    //Dictionary<int, IMatrix2D>  ks = new Dictionary<int, IMatrix2D>(model.SubdomainsDictionary.Count);
        //    Dictionary<int, double[][]> KfpDqSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
        //    Dictionary<int, double[][]> KppDqVectorsSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);

        //    int procs = VectorExtensions.AffinityCount;
        //    var k = model.SubdomainsDictionary.Keys.Select(x => x).ToArray<int>();
        //    var internalKs = new Dictionary<int, IMatrix2D>[procs];
        //    var internalKfpDq = new Dictionary<int, double[][]>[procs];
        //    var internalKppDq = new Dictionary<int, double[][]>[procs];

        //    //foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)" tou opoiu ta stoixeia ginontai access kai me ID
        //    //{
        //    //    KfpDqSubdomains.Add(subdomain.ID, SubdomainCalculations.CalculateKfreeprescribedDqMultiplications(subdomain, elementProvider, scaleTransitions, boundaryNodes));
        //    //}

        //    Parallel.ForEach(k.PartitionLimits(procs), limit =>
        //    {
        //        if (limit.Item3 - limit.Item2 > 0)
        //        {
        //            internalKs[limit.Item1] = new Dictionary<int, IMatrix2D>(limit.Item3 - limit.Item2);
        //            internalKfpDq[limit.Item1] = new Dictionary<int, double[][]>(limit.Item3 - limit.Item2);
        //            internalKppDq[limit.Item1] = new Dictionary<int, double[][]>(limit.Item3 - limit.Item2);
        //            for (int i = limit.Item2; i < limit.Item3; i++)
        //            {
        //                //SkylineMatrix2D mat = GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(model.SubdomainsDictionary[k[i]], s);

        //                // prosthiki print
        //                //ekteleseis_counter += 1;
        //                //string counter_data = ekteleseis_counter.ToString();
        //                //string file_no = (k[i]).ToString();
        //                //string print_path = string.Format(print_path_gen, file_no, counter_data);
        //                //mat.WriteToFile(print_path);
        //                (double[][]KfpDqVectors,  double[][]KppDqVectors, SkylineMatrix2D K)= UpdateSubdomainKffAndCalculateKfpDqAndKppDqp(model.SubdomainsDictionary[k[i]], elementProvider,scaleTransitions,
        //                boundaryNodes, boundaryElements[k[i]]);
        //                internalKs[limit.Item1].Add(k[i], K);
        //                internalKfpDq[limit.Item1].Add(k[i], KfpDqVectors); 
        //                internalKppDq[limit.Item1].Add(k[i], KppDqVectors);  
        //            }
        //        }
        //        else
        //        {
        //            internalKs[limit.Item1] = new Dictionary<int, IMatrix2D>();
        //            internalKfpDq[limit.Item1] = new Dictionary<int, double[][]>();
        //            internalKppDq[limit.Item1] = new Dictionary<int, double[][]>();
        //        }
        //    });

        //    for (int i = 0; i < procs; i++)
        //    {
        //        foreach (int key in internalKs[i].Keys)
        //        {
        //            //ks.Add(key, internalKs[i][key]);
        //            linearSystems[key].Matrix = internalKs[i][key];
        //            KfpDqSubdomains.Add(key, internalKfpDq[i][key]);
        //            KppDqVectorsSubdomains.Add(key, internalKppDq[i][key]);
        //        }
        //    }


        //    //(Dictionary<int, double[][]> KfpDqSubdomains,Dictionary<int, double[][]> KppDqVectorsSubdomains)

        //    return (KfpDqSubdomains, KppDqVectorsSubdomains);
        //}

        private double[][] KfpDqVectors;
        private double[][] KppDqVectors;
        ISubdomainFreeDofOrdering dofOrdering;
        DofTable FreeDofs;
        //v2.1 Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary;
        IScaleTransitions_v2 scaleTransitions;
        Dictionary<int, Dictionary<int, Element_v2>> boundaryElements;
        Dictionary<int, Node_v2> boundaryNodes;
        Dictionary<int, int> boundaryNodesOrder;
        int currentSubdomainID;

        public (Dictionary<int, double[][]>, Dictionary<int, double[][]>) UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultipleObje_v2(Model_v2 model, IElementMatrixProvider_v2 elementProvider, IScaleTransitions_v2 scaleTransitions,
            Dictionary<int, Node_v2> boundaryNodes, Dictionary<int, Dictionary<int, Element_v2>> boundaryElements,ISolver_v2 solver)
        {
            IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems = solver.LinearSystems; //v2.3

            Dictionary<int, double[][]> KfpDqSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            Dictionary<int, double[][]> KppDqVectorsSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            this.boundaryElements = boundaryElements;
            this.boundaryNodes = boundaryNodes;
            this.scaleTransitions = scaleTransitions;

            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                dofOrdering = subdomain.DofOrdering; //_v2.1
                FreeDofs = subdomain.DofOrdering.FreeDofs;//_v2.1 nodalDOFsDictionary = subdomain.NodalDOFsDictionary;
                currentSubdomainID = subdomain.ID;

                #region Create KfpDq and KppDq vectors 
                KfpDqVectors = new double[scaleTransitions.MacroscaleVariableDimension()][];
                for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
                {
                    KfpDqVectors[j1] = new double[dofOrdering.NumFreeDofs]; //v2.2 subdomain.TotalDOFs]; 
                }

                KppDqVectors = new double[scaleTransitions.MacroscaleVariableDimension()][];
                boundaryNodesOrder = SubdomainCalculations_v2.GetNodesOrderInDictionary(boundaryNodes);
                for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
                {
                    KppDqVectors[j1] = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)
                }
                #endregion

                var StiffnessProvider = new StiffnessProviderSimu_v2(this);

                var subdomainK = solver.BuildGlobalMatrix(subdomain, StiffnessProvider);
                //v2.4 var subdomainK= GlobalMatrixAssemblerSkyline.CalculateFreeFreeGlobalMatrix(subdomain, StiffnessProvider);                

                linearSystems[subdomain.ID].Matrix=subdomainK;
                //v2.5 linearSystems[subdomain.ID].Matrix = subdomainK;

                KfpDqSubdomains.Add(subdomain.ID, KfpDqVectors);
                KppDqVectorsSubdomains.Add(subdomain.ID, KppDqVectors);                
            }

            return (KfpDqSubdomains, KppDqVectorsSubdomains);        
        }

        public void UpdateVectors_v2(IElement_v2 element, IMatrix ElementK)
        {
            if (boundaryElements[currentSubdomainID].ContainsKey(element.ID))//COPIED From UpdateSubdomainKffAndCalculateKfpDqAndKppDqp (prosoxh boundary elements Dictionary diathetoun kai to model kai to subdomain kai einai diaforetika edw exei diorthwthei
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
                        bool isFree = FreeDofs.TryGetValue(matrixAssemblyNodes[i], elementDOFTypes[i][dofTypeRowToNumber],
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
                                        KfpDqVectors[j2][dofRow] += contribution[j2]; // TODO diorthothike
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
                                        KppDqVectors[j1][dofrow_p] += contribution[j1];
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
