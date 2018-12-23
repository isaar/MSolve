using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: The F = Ff - Kfc*Fc should not be done in the solver. The solver should only operate on the final linear systems.
//      It could be done here or in the analyzer.
//TODO: Optimizations should be available when building matrices with the same sparsity pattern. These optimizations should work
//      for all assemblers, not only skyline.
//TODO: remove casts and the logic for specific elements
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineAssembler : IGlobalMatrixAssembler<SkylineMatrix>
    {
        private const string name = "SkylineAssembler"; // for error messages
        //private readonly LinearSystem_v2<SkylineMatrix, Vector> linearSystem;

        //public SkylineAssembler(IReadOnlyList<LinearSystem_v2<SkylineMatrix, LegacyVector>> linearSystems)
        //{
        //    if (linearSystems.Count != 1) throw new InvalidMatrixFormatException(
        //        name + " can be used if there is only 1 subdomain.");
        //    this.linearSystem = linearSystems[0];
        //}

        //TODO: This is for the case when we also number constrained dofs globally.
        //public (SkylineMatrix Kff, DokRowMajor Kfc) BuildGlobalMatrices(IEnumerable<IElement> elements,
        //    AllDofOrderer dofOrderer, IElementMatrixProvider matrixProvider)
        //{
        //    int numConstrainedDofs = dofOrderer.NumConstrainedDofs;
        //    int numFreeDofs = dofOrderer.NumFreeDofs;
        //    SkylineBuilder Kff = FindSkylineColumnHeights(elements, numFreeDofs, dofOrderer.FreeDofs);
        //    var Kfc = DokRowMajor.CreateEmpty(numFreeDofs, numConstrainedDofs);

        //    foreach (IElement elementWrapper in elements)
        //    {
        //        ContinuumElement2D element = (ContinuumElement2D)(elementWrapper.IElementType);
        //        // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
        //        (IReadOnlyDictionary<int, int> mapStandard, IReadOnlyDictionary<int, int> mapConstrained) =
        //            dofOrderer.MapDofsElementToGlobal(element);
        //        Matrix k = Conversions.MatrixOldToNew(matrixProvider.Matrix(elementWrapper));
        //        Kff.AddSubmatrixSymmetric(k, mapStandard);
        //        Kfc.AddSubmatrix(k, mapStandard, mapConstrained);
        //    }

        //    //TODO: perhaps I should filter the matrices in the concrete class before returning (e.g. dok.Build())
        //    return (Kff.BuildSkylineMatrix(), Kfc);
        //}

        public SkylineMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements, 
            IElementMatrixProvider matrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            SkylineBuilder Kff = FindSkylineColumnHeights(elements, numFreeDofs, dofOrdering.FreeDofs);

            foreach (IElement element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                //IReadOnlyDictionary<int, int> elementToGlobalDofs = dofOrdering.MapFreeDofsElementToSubdomain(element);
                Matrix k = matrixProvider.Matrix(element).LegacyToNewMatrix();
                Kff.AddSubmatrixSymmetric(k, elementDofIndices, subdomainDofIndices);
            }

            return Kff.BuildSkylineMatrix();
        }

        public SkylineMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider matrixProvider) //TODO: remove this
        {
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;
            var K = new Numerical.LinearAlgebra.SkylineMatrix2D(CalculateRowIndex(subdomain));
            foreach (IElement element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.IElementType is FEM.Interfaces.IEmbeddedElement;
                IMatrix2D ElementK = matrixProvider.Matrix(element);

                var elementDOFTypes = element.IElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    foreach (DOFType dofTypeRow in elementDOFTypes[i])
                    {
                        int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (dofRow != -1)
                        {
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                {
                                    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                    if (dofColumn != -1)
                                    {
                                        int height = dofRow - dofColumn;
                                        if (height >= 0)
                                            K.Data[K.RowIndex[dofRow] + height] += ElementK[iElementMatrixRow, iElementMatrixColumn];
                                    }
                                    iElementMatrixColumn++;
                                }
                            }
                        }
                        iElementMatrixRow++;
                    }
                }
            }

            return SkylineMatrix.CreateFromArrays(K.Columns, K.Data, K.RowIndex, false);
        }

        //TODO: If one element engages some dofs (of a node) and another engages other dofs, the ones not in the intersection 
        // are not dependent from the rest. This method assumes dependency for all dofs of the same node. This is a rare occasion 
        // though.
        private static SkylineBuilder FindSkylineColumnHeights(IEnumerable<IElement> elements,
            int numFreeDofs, DofTable freeDofs)
        {
            int[] colHeights = new int[numFreeDofs]; //only entries above the diagonal count towards the column height
            foreach (IElement element in elements)
            {
                //TODO: perhaps the 2 outer loops could be done at once to avoid a lot of dof indexing. Could I update minDof
                //      and colHeights[] at once?

                IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);

                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0.0 entries in the element stiffness matrix.
                int minDof = Int32.MaxValue;
                foreach (var node in elementNodes)
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node)) minDof = Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (var node in elementNodes)
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return SkylineBuilder.Create(numFreeDofs, colHeights);
        }

        //TODO: If one element engages some dofs (of a node) and another engages other dofs, the ones not in the intersection 
        // are not dependent from the rest. This method assumes dependency for all dofs of the same node. This is a rare occasion 
        // though.
        private static int[] CalculateRowIndex(ISubdomain subdomain) //TODO: remove this
        {
            int order = subdomain.TotalDOFs;
            int[] rowHeights = new int[order];
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

            foreach (IElement element in subdomain.ElementsDictionary.Values)
            {
                int minDOF = Int32.MaxValue;
                foreach (INode node in element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.IElementType is FEM.Interfaces.IEmbeddedElement))
                        continue;
                    foreach (int dof in nodalDOFsDictionary[node.ID].Values)
                        if (dof != -1) minDOF = Math.Min(dof, minDOF);
                }
                //foreach (Node node in element.NodesDictionary.Values.Where(e => e.EmbeddedInElement == null))
                foreach (INode node in element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.IElementType is FEM.Interfaces.IEmbeddedElement))
                        continue;
                    foreach (int dof in nodalDOFsDictionary[node.ID].Values)
                        if (dof != -1) rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                }
            }

            int[] rowIndex = new int[subdomain.TotalDOFs + 1];
            rowIndex[0] = 0;
            rowIndex[1] = 1;
            for (int i = 1; i < subdomain.TotalDOFs; i++) rowIndex[i + 1] = rowIndex[i] + rowHeights[i] + 1;

            return rowIndex;
        }
    }
}
