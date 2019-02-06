using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System.Collections;
using System.Linq;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;

namespace ISAAR.MSolve.Analyzers.SupportiveClasses
{
    public static class SubdomainEqForcesCalculations
    {
        public static Vector CalculateKfreeprescribedUpMultiplicationForSubdRHSContribution(Subdomain subdomain,
            Dictionary<int, Node> boundaryNodes,Dictionary<int,Dictionary<DOFType,double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            int nIncrement,int totalIncrements)
        {
            ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider(); //TODOMS: where the provider should be defined
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

            double[] Kfp_Ustep  = new double[subdomain.TotalDOFs]; // h allliws subdomain.Forces.GetLength(0)
            
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);

            foreach (Element element in subdomain.ElementsDictionary.Values) // TODOGerasimos edw mporei na xrhsimopoihthei to dictionary twn eleement pou exoun fp nodes
            {
                var isEmbeddedElement = element.ElementType is IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    foreach (DOFType dofTypeRow in elementDOFTypes[i])
                    {
                        int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (dofRow != -1) // TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
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
                                    //foreach(KeyValuePair< DOFType,double>  prescribedDOFtype in totalBoundaryDisplacements[nodeColumn.ID])
                                    //{
                                    //    prescribedDOFtype.Key
                                    //}

                                    double[] element_Kfp_triplette = new double[nodalDofsNumber]; //nodalDofsNumber: giati oxi scaleTransitions.PrescribedDofsPerNode()? Dioti tou ta pairname ola(triplette) kai dialegei to 
                                    for (int j1 = 0; j1 < nodalDofsNumber; j1++)                  //scaleTransitions.MicroToMacroTransition ti tha xrhsimopoihsei apo afta analoga pws einai implemented
                                    {
                                        element_Kfp_triplette[j1] = ElementK[iElementMatrixRow, iElementMatrixColumn + j1];
                                    }


                                    Dictionary<DOFType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
                                    Dictionary<DOFType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
                                    double[] uStep_values_orZero_for_free = new double[nodalDofsNumber];

                                    int positionOfDof =0;
                                    foreach( DOFType doftype1 in elementDOFTypes[j] )
                                    {
                                        if (nodalConvergedDisplacements.ContainsKey(doftype1))
                                        {
                                            uStep_values_orZero_for_free[positionOfDof] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
                                            // TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
                                        }
                                        positionOfDof += 1;
                                    }

                                    double contribution = 0;
                                    for (int j2=0; j2<nodalDofsNumber;j2++)
                                    {
                                        contribution += element_Kfp_triplette[j2] * uStep_values_orZero_for_free[j2];
                                    }

                                    Kfp_Ustep[dofRow] += contribution;

                                    

                                }
                                iElementMatrixColumn += nodalDofsNumber;

                            }
                        }
                        iElementMatrixRow++;
                    }
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return new Vector(Kfp_Ustep);
        }
    }
}
