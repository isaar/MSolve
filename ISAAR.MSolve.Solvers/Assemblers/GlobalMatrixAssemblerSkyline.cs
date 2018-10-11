using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using IEmbeddedElement = ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement; // TODO: Probably solvers should not have dependencies from FEM, IGA, etc. Find a way to hide the embedding.

namespace ISAAR.MSolve.Solvers.Assemblers
{
    public static class GlobalMatrixAssemblerSkyline
    {
        private static int[] CalculateRowIndex(ISubdomain subdomain, Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary)
        {
            int[] rowHeights = new int[subdomain.TotalDOFs];
            foreach (IElement element in subdomain.ΙElementsDictionary.Values)
            {
                int minDOF = Int32.MaxValue;
                foreach (INode node in element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.IElementType is IEmbeddedElement))
                        continue;
                    foreach (int dof in nodalDOFsDictionary[node.ID].Values)
                        if (dof != -1) minDOF = Math.Min(dof, minDOF);
                }
                //foreach (Node node in element.NodesDictionary.Values.Where(e => e.EmbeddedInElement == null))
                foreach (INode node in element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.IElementType is IEmbeddedElement))
                        continue;
                    foreach (int dof in nodalDOFsDictionary[node.ID].Values)
                        if (dof != -1) rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                }
            }

            int[] rowIndex = new int[subdomain.TotalDOFs + 1];
            rowIndex[0] = 0;
            rowIndex[1] = 1;
            for (int i = 1; i < subdomain.TotalDOFs; i++)
                rowIndex[i + 1] = rowIndex[i] + rowHeights[i] + 1;
            return rowIndex;
        }
		
        public static SkylineMatrix2D CalculateGlobalMatrix(ISubdomain subdomain, Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary, IElementMatrixProvider elementProvider)
        {
            // TODO: should encapsulate DOF logic into a separate entity that will manage things if embedded or not (should return element matrix and globaldofs correspondence list
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            SkylineMatrix2D K = new SkylineMatrix2D(GlobalMatrixAssemblerSkyline.CalculateRowIndex(subdomain, nodalDOFsDictionary));
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (IElement element in subdomain.ΙElementsDictionary.Values)
            {
                var isEmbeddedElement = element.IElementType is IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
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
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;
            return K;
        }

        public static SkylineMatrix2D CalculateGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementProvider)
        {
            return CalculateGlobalMatrix(subdomain, subdomain.NodalDOFsDictionary, elementProvider);
        }
    }
}
