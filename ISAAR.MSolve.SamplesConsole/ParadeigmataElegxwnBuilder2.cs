//using ISAAR.MSolve.PreProcessor;
//using ISAAR.MSolve.PreProcessor.Elements;
//using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
//using ISAAR.MSolve.PreProcessor.Embedding;
// compa
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;//using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
//using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.SamplesConsole
{
    class ParadeigmataElegxwnBuilder2
    {
        public static void ShellPlateBuilder(Model model, double load_value)
        {
            // proelefsi: branch master idio onoma ParadeigmataElegxwnBuilder2.ShellPlateBuilder(Model model, double load_value)
            ShellElasticMaterial material1 = new ShellElasticMaterial()
            {
                YoungModulus = 135300,
                PoissonRatio = 0.3,
                ShearCorrectionCoefficientK = 5 / 6,
            };

            double[,] nodeData = new double[,] {{10.000000,10.000000,0.000000},
            {7.500000,10.000000,0.000000},
            {5.000000,10.000000,0.000000},
            {2.500000,10.000000,0.000000},
            {0.000000,10.000000,0.000000},
            {10.000000,7.500000,0.000000},
            {7.500000,7.500000,0.000000},
            {5.000000,7.500000,0.000000},
            {2.500000,7.500000,0.000000},
            {0.000000,7.500000,0.000000},
            {10.000000,5.000000,0.000000},
            {7.500000,5.000000,0.000000},
            {5.000000,5.000000,0.000000},
            {2.500000,5.000000,0.000000},
            {0.000000,5.000000,0.000000},
            {10.000000,2.500000,0.000000},
            {7.500000,2.500000,0.000000},
            {5.000000,2.500000,0.000000},
            {2.500000,2.500000,0.000000},
            {0.000000,2.500000,0.000000},
            {10.000000,0.000000,0.000000},
            {7.500000,0.000000,0.000000},
            {5.000000,0.000000,0.000000},
            {2.500000,0.000000,0.000000},
            {0.000000,0.000000,0.000000},
            {8.750000,10.000000,0.000000},
            {6.250000,10.000000,0.000000},
            {3.750000,10.000000,0.000000},
            {1.250000,10.000000,0.000000},
            {8.750000,7.500000,0.000000},
            {6.250000,7.500000,0.000000},
            {3.750000,7.500000,0.000000},
            {1.250000,7.500000,0.000000},
            {8.750000,5.000000,0.000000},
            {6.250000,5.000000,0.000000},
            {3.750000,5.000000,0.000000},
            {1.250000,5.000000,0.000000},
            {8.750000,2.500000,0.000000},
            {6.250000,2.500000,0.000000},
            {3.750000,2.500000,0.000000},
            {1.250000,2.500000,0.000000},
            {8.750000,0.000000,0.000000},
            {6.250000,0.000000,0.000000},
            {3.750000,0.000000,0.000000},
            {1.250000,0.000000,0.000000},
            {10.000000,8.750000,0.000000},
            {10.000000,6.250000,0.000000},
            {10.000000,3.750000,0.000000},
            {10.000000,1.250000,0.000000},
            {7.500000,8.750000,0.000000},
            {7.500000,6.250000,0.000000},
            {7.500000,3.750000,0.000000},
            {7.500000,1.250000,0.000000},
            {5.000000,8.750000,0.000000},
            {5.000000,6.250000,0.000000},
            {5.000000,3.750000,0.000000},
            {5.000000,1.250000,0.000000},
            {2.500000,8.750000,0.000000},
            {2.500000,6.250000,0.000000},
            {2.500000,3.750000,0.000000},
            {2.500000,1.250000,0.000000},
            {0.000000,8.750000,0.000000},
            {0.000000,6.250000,0.000000},
            {0.000000,3.750000,0.000000},
            {0.000000,1.250000,0.000000}, };

            int[,] elementData = new int[,] { {1,1,2,7,6,26,50,30,46},
            {2,2,3,8,7,27,54,31,50},
            {3,3,4,9,8,28,58,32,54},
            {4,4,5,10,9,29,62,33,58},
            {5,6,7,12,11,30,51,34,47},
            {6,7,8,13,12,31,55,35,51},
            {7,8,9,14,13,32,59,36,55},
            {8,9,10,15,14,33,63,37,59},
            {9,11,12,17,16,34,52,38,48},
            {10,12,13,18,17,35,56,39,52},
            {11,13,14,19,18,36,60,40,56},
            {12,14,15,20,19,37,64,41,60},
            {13,16,17,22,21,38,53,42,49},
            {14,17,18,23,22,39,57,43,53},
            {15,18,19,24,23,40,61,44,57},
            {16,19,20,25,24,41,65,45,61},};

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });
            }

            // orismos elements 
            Element e1;
            int subdomainID = 1;
            double tk_shell_plate = 0.5;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Shell8dispCopyGetRAM_1(material1, 3, 3, 2)//ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
                    {
                        //oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
                        oVn_i = new double[][] { new double[] { 0,0,1 },
                                                 new double[] { 0,0,1 },
                                                 new double[] { 0,0,1 },
                                                 new double[] { 0,0,1 },
                                                 new double[] { 0,0,1 },
                                                 new double[] { 0,0,1 },
                                                 new double[] { 0,0,1 },
                                                 new double[] { 0,0,1 },},
                        tk = new double[] { tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate, tk_shell_plate },
                    }
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            }

            // constraint paaktwsh gurw gurw plevres
            int pointID;
            int[] cnstrnd = new int[] { 21, 22, 23, 24, 25, 26, 27, 28, 29, 1, 2, 3, 4, 5, 42, 43, 44, 45, 46, 47, 48, 49, 6, 11, 16, 10, 15, 20, 62, 63, 64, 65 };
            for (int k = 0; k < cnstrnd.GetLength(0); k++)
            {
                pointID = cnstrnd[k];
                model.NodesDictionary[pointID].Constraints.Add(DOFType.X);
                model.NodesDictionary[pointID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[pointID].Constraints.Add(DOFType.Z);
                model.NodesDictionary[pointID].Constraints.Add(DOFType.RotX);
                model.NodesDictionary[pointID].Constraints.Add(DOFType.RotY);
                
            }

            // fortish korufhs
            Load load1;
            
            load1 = new Load()
            {
                    Node = model.NodesDictionary[13],
                    DOF = DOFType.Z,
                    Amount = 1 * load_value
            };
            model.Loads.Add(load1);
            
        }
    }
}
