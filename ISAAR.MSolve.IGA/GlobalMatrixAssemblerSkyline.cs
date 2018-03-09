using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA
{
    public static class GlobalMatrixAssemblerSkyline
    {
        private static int[] CalculateRowIndex(Patch patch, Dictionary<int, Dictionary<DOFType,int>> controlPointDOFsDictionary)
        {
            int[] rowHeights = new int[patch.TotalDOFs];
            foreach (Element element in patch.ElementsDictionary.Values)
            {
                int minDOF = Int32.MaxValue;
                foreach (ControlPoint controlPoint in element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    foreach (int dof in controlPointDOFsDictionary[controlPoint.ID].Values)
                        if (dof != -1)
                            minDOF = Math.Min(dof, minDOF);
                }

                foreach (ControlPoint controlPoint in element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element))
                {
                    foreach (int dof in controlPointDOFsDictionary[controlPoint.ID].Values)
                        if (dof != -1)
                            rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                }
            }

            int[] rowIndex = new int[patch.TotalDOFs + 1];
            rowIndex[0] = 0;
            rowIndex[1] = 1;
            for (int i = 1; i < patch.TotalDOFs; i++)
                rowIndex[i + 1] = rowIndex[i] + rowHeights[i] + 1;
            return rowIndex;
        }

        public static SkylineMatrix2D CalculateGlobalMatrix(Patch patch, Dictionary<int, Dictionary<DOFType, int>> controlPointDOFsDictionary, IElementMatrixProvider elementProvider)
        {
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            SkylineMatrix2D K = new SkylineMatrix2D(GlobalMatrixAssemblerSkyline.CalculateRowIndex(patch, controlPointDOFsDictionary));
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element element in patch.ElementsDictionary.Values)
            {
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    ControlPoint controPointRow = matrixAssemblyNodes[i];
                    foreach (DOFType dofTypeRow in elementDOFTypes[i])
                    {
                        /////// Ti ypologizei to dofRow
                        int dofRow = controlPointDOFsDictionary.ContainsKey(controPointRow.ID) == false ? -1 :controlPointDOFsDictionary[controPointRow.ID][dofTypeRow];
                        ;
                        if (dofRow != -1)
                        {
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                ControlPoint controlPointColumn = matrixAssemblyNodes[j];
                                foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                {
                                    int dofColumn = controlPointDOFsDictionary.ContainsKey(controlPointColumn.ID) == false ? -1 : controlPointDOFsDictionary[controlPointColumn.ID][dofTypeColumn];
                                    ;
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


        public static SkylineMatrix2D CalculateGlobalMatrix(Patch patch, IElementMatrixProvider elementProvider)
        {
            return CalculateGlobalMatrix(patch, patch.ControlPointDOFsDictionary, elementProvider);
        }
    }
}
