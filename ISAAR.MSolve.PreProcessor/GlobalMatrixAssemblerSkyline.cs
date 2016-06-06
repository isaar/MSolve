using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.PreProcessor.Interfaces;
using System.IO;

namespace ISAAR.MSolve.PreProcessor
{
    public static class GlobalMatrixAssemblerSkyline
    {
        private static int[] CalculateRowIndex(Subdomain subdomain, Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary)
        {
            int[] rowHeights = new int[subdomain.TotalDOFs];
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                int minDOF = Int32.MaxValue;
                //foreach (Node node in element.NodesDictionary.Values.Where(e => e.EmbeddedInElement == null))
                foreach (Node node in element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.ElementType is IEmbeddedElement))
                        continue;
                    foreach (int dof in nodalDOFsDictionary[node.ID].Values)
                        if (dof != -1) minDOF = Math.Min(dof, minDOF);
                }
                //foreach (Node node in element.NodesDictionary.Values.Where(e => e.EmbeddedInElement == null))
                foreach (Node node in element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    if ((nodalDOFsDictionary.ContainsKey(node.ID) == false) && (element.ElementType is IEmbeddedElement))
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

        //private static int[] CalculateRowIndex(Subdomain subdomain)
        //{
        //    return CalculateRowIndex(subdomain, subdomain.NodalDOFsDictionary);
        //}

        public static SkylineMatrix2D<double> CalculateGlobalMatrix(Subdomain subdomain, Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary, IElementMatrixProvider elementProvider)
        {
            // TODO: should encapsulate DOF logic into a separate entity that will manage things if embedded or not (should return element matrix and globaldofs correspondence list
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            SkylineMatrix2D<double> K = new SkylineMatrix2D<double>(GlobalMatrixAssemblerSkyline.CalculateRowIndex(subdomain, nodalDOFsDictionary));
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D<double> ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    Node nodeRow = matrixAssemblyNodes[i];
                    foreach (DOFType dofTypeRow in elementDOFTypes[i])
                    {
                        int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (dofRow != -1)
                        {
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                Node nodeColumn = matrixAssemblyNodes[j];
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

        public static SkylineMatrix2D<double> CalculateGlobalMatrix(Subdomain subdomain, IElementMatrixProvider elementProvider)
        {
            return CalculateGlobalMatrix(subdomain, subdomain.NodalDOFsDictionary, elementProvider);
        }
    }
}
