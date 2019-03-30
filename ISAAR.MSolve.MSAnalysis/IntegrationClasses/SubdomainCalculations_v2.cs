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
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Matrices;
//using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM
{
    /// <summary>
    /// Supportive class  that implements nesessary integration methods associated with FE2 multiscale analysis 
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class SubdomainCalculations_v2
    {
        
        public static Dictionary<int,int> GetNodesOrderInDictionary(Dictionary<int, Node_v2> boundaryNodes)
        {
            Dictionary<int, int> boundaryNodesOrder = new Dictionary<int, int>();
            int order = 1;
            foreach (Node_v2 boundaryNode in boundaryNodes.Values)
            {
                boundaryNodesOrder.Add(boundaryNode.ID, order);
                order += 1;
            }
            return boundaryNodesOrder;

        }

        
        #region v2 methods
        public static double[][] SubtractConsecutiveVectors_v2(double[][] KppDqVectors, double[][] f3_vectors)
        {
            double[][] f4_vectors = new double[KppDqVectors.GetLength(0)][];
            for (int j1 = 0; j1 < f4_vectors.GetLength(0); j1++)
            {
                f4_vectors[j1] = new double[KppDqVectors[0].GetLength(0)];
                for (int j2 = 0; j2 < f4_vectors[0].GetLength(0); j2++)
                {
                    f4_vectors[j1][j2] = KppDqVectors[j1][j2] - f3_vectors[j1][j2];
                }
            }

            return f4_vectors;

        }

        public static double[,] CalculateDqCondDq_v2(double[][] f4_vectors, IScaleTransitions_v2 scaleTransitions, Dictionary<int, Node_v2> boundaryNodes)
        {
            double[,] DqCondDq = new double[scaleTransitions.MacroscaleVariableDimension(), scaleTransitions.MacroscaleVariableDimension()];

            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);

            foreach (Node_v2 boundaryNode in boundaryNodes.Values)
            {
                for (int i1 = 0; i1 < f4_vectors.GetLength(0); i1++)
                {
                    double[] f4DataTriplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                    for (int i2 = 0; i2 < scaleTransitions.PrescribedDofsPerNode(); i2++)
                    {
                        f4DataTriplette[i2] = f4_vectors[i1][scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[boundaryNode.ID] - 1) + i2];
                    }
                    double[] contribution = scaleTransitions.MicroToMacroTransition(boundaryNode, f4DataTriplette);
                    for (int i3 = 0; i3 < scaleTransitions.MacroscaleVariableDimension(); i3++)
                    {
                        DqCondDq[i3, i1] += contribution[i3];
                    }
                }

            }
            return DqCondDq;
        }

        public static double[] CalculateDqFpp_v2(double[] FppReactionVector, IScaleTransitions_v2 scaleTransitions, Dictionary<int, Node_v2> boundaryNodes)
        {
            double[] DqFpp = new double[scaleTransitions.MacroscaleVariableDimension()];

            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);

            foreach (Node_v2 boundaryNode in boundaryNodes.Values)
            {

                double[] FppDataTriplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                for (int i2 = 0; i2 < scaleTransitions.PrescribedDofsPerNode(); i2++)
                {
                    FppDataTriplette[i2] = FppReactionVector[scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[boundaryNode.ID] - 1) + i2];
                }
                double[] contribution = scaleTransitions.MicroToMacroTransition(boundaryNode, FppDataTriplette);
                for (int i3 = 0; i3 < scaleTransitions.MacroscaleVariableDimension(); i3++)
                {
                    DqFpp[i3] += contribution[i3];
                }


            }

            return DqFpp;

        }

        public static double[] CalculateFppReactionsVector_v2(Subdomain_v2 subdomain, IElementMatrixProvider_v2 elementProvider,
            IScaleTransitions_v2 scaleTransitions, Dictionary<int, Node_v2> boundaryNodes, IVectorView solution, IVectorView dSolution,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {
            //TODOGerasimos: 1) Subdomain2 einai h upo kataskevh subdomain.cs ths Marias gia na mporoume na anaferthoume sthn methodo ths CalculateElementNodalDisplacements(..,..). 
            // Otan parei telikh morfh tha taftizetai me thn Subdomain.cs
            // 2)IVector solution, IVector dSolution EINAI AFTA ME TA OPOIA kaloume thn GetRHSFromSolution sthn 213 tou NRNLAnalyzer
            double[] FppReactionVector;
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);
            FppReactionVector = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element_v2 element in subdomain.Elements)
            {
                var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
                var elStart = DateTime.Now;
                //IMatrix2D ElementK = elementProvider.Matrix(element);
                var localSolution = subdomain.GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, solution);
                subdomain.ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
                double[] localdSolution = subdomain.GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, dSolution);
                double[] f = element.ElementType.CalculateForces(element, localSolution, localdSolution);

                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DofEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodes.ContainsKey(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            FppReactionVector[dofrow_p] += f[iElementMatrixRow + i1];

                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return FppReactionVector;

        }

        public static double[][] CalculateKpfKffinverseKfpDq_v2(double[][] f2_vectors, Subdomain_v2 subdomain, IElementMatrixProvider_v2 elementProvider, IScaleTransitions_v2 scaleTransitions, Dictionary<int, Node_v2> boundaryNodes)
        {
            var dofOrdering = subdomain.FreeDofOrdering; //_v2.1
            var FreeDofs = subdomain.FreeDofOrdering.FreeDofs;//_v2.1Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

            double[][] f3_vectors = new double[f2_vectors.GetLength(0)][];
            for (int i1 = 0; i1 < f2_vectors.GetLength(0); i1++)
            {
                f3_vectors[i1] = new double[scaleTransitions.PrescribedDofsPerNode() * boundaryNodes.Count];
            }
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element_v2 element in subdomain.Elements) //_v2.3 ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DofEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
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
                                int dofTypeColumnToNumber = -1;
                                foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                {
                                    dofTypeColumnToNumber++;
                                    bool isFree = FreeDofs.TryGetValue(matrixAssemblyNodes[j], elementDOFTypes[j][dofTypeColumnToNumber],
                                    out int dofColumn); // v2.4 int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                    if (isFree)// TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                                    {                   // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok

                                        for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
                                        {
                                            f3_vectors[i2][dofrow_p] += ElementK[iElementMatrixRow + i1, iElementMatrixColumn] * f2_vectors[i2][dofColumn];
                                        }

                                    }
                                    iElementMatrixColumn++;
                                }
                            }

                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return f3_vectors;

        }
        #endregion
    }
}
