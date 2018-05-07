using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Tests.Supportive_Classes
{
    public static class HexaCantileverBeamForNewtonRaphson
    {
        public static void Build(Model model, double load_value)
        {
            var material1 = new ElasticMaterial3D()
            {
                YoungModulus = 1353000,
                PoissonRatio = 0.3,
            };

            var nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
            {0.250000,-0.250000,-1.000000},
            {-0.250000,0.250000,-1.000000},
            {0.250000,0.250000,-1.000000},
            {-0.250000,-0.250000,-0.500000},
            {0.250000,-0.250000,-0.500000},
            {-0.250000,0.250000,-0.500000},
            {0.250000,0.250000,-0.500000},
            {-0.250000,-0.250000,0.000000},
            {0.250000,-0.250000,0.000000},
            {-0.250000,0.250000,0.000000},
            {0.250000,0.250000,0.000000},
            {-0.250000,-0.250000,0.500000},
            {0.250000,-0.250000,0.500000},
            {-0.250000,0.250000,0.500000},
            {0.250000,0.250000,0.500000},
            {-0.250000,-0.250000,1.000000},
            {0.250000,-0.250000,1.000000},
            {-0.250000,0.250000,1.000000},
            {0.250000,0.250000,1.000000}};

            var elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
            {2,12,11,9,10,8,7,5,6},
            {3,16,15,13,14,12,11,9,10},
            {4,20,19,17,18,16,15,13,14}, };

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

            }

            // orismos elements 
            Element e1;
            var subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element
                {
                    ID = nElement + 1,
                    //ElementType = new Hexa8NL(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                    //ElementType = new Hexa8(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };

                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            }

            // constraint vashh opou z=-1
            for (int k = 1; k < 5; k++)
            {
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Z });
            }

            // fortish korufhs
            Load load1;
            for (int k = 17; k < 21; k++)
            {
                load1 = new Load
                {
                    Node = model.NodesDictionary[k],
                    DOF = DOFType.X,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }
    }
}
