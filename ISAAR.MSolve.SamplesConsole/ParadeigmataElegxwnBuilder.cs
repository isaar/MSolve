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
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
//using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.SamplesConsole
{
    class ParadeigmataElegxwnBuilder
    {
        //public static void Example_cohesive_hexa_mixed(Model model)
        //{
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8node(model); 
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeConstraintsMixed(model);
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeLoadsMIxed(model, 1);
        //}

        //public static void Example_cohesive_hexa_orthi_elastic(Model model)
        //{
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8node(model); // me 1353000
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeConstraintsOrthiElastic(model);
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeLoadsElasticOrthi(model, 3.47783);
        //}

        //public static void Example_cohesive_hexa_orthi_constr_anw_benc1(Model model)
        //{
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8node(model); // me 135300
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeConstraintsBenc1(model);
        //    ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeLoadsBenc1(model, 8);// gia elastiko klado 
        //    //ParadeigmataElegxwnBuilder.Example2Hexa8NL1Cohesive8nodeLoadsBenc1(model, 13.75); // gia metelastiko
        //}

        #region to comment out temporarily

        //public static void ExampleHexaCantilever(Model model)
        //{

        //}

        //public static void HexaCantileverBuilder(Model model, double load_value)
        //{
        //    ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };

        //    double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
        //    {0.250000,-0.250000,-1.000000},
        //    {-0.250000,0.250000,-1.000000},
        //    {0.250000,0.250000,-1.000000},
        //    {-0.250000,-0.250000,-0.500000},
        //    {0.250000,-0.250000,-0.500000},
        //    {-0.250000,0.250000,-0.500000},
        //    {0.250000,0.250000,-0.500000},
        //    {-0.250000,-0.250000,0.000000},
        //    {0.250000,-0.250000,0.000000},
        //    {-0.250000,0.250000,0.000000},
        //    {0.250000,0.250000,0.000000},
        //    {-0.250000,-0.250000,0.500000},
        //    {0.250000,-0.250000,0.500000},
        //    {-0.250000,0.250000,0.500000},
        //    {0.250000,0.250000,0.500000},
        //    {-0.250000,-0.250000,1.000000},
        //    {0.250000,-0.250000,1.000000},
        //    {-0.250000,0.250000,1.000000},
        //    {0.250000,0.250000,1.000000}};

        //    int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
        //    {2,12,11,9,10,8,7,5,6},
        //    {3,16,15,13,14,12,11,9,10},
        //    {4,20,19,17,18,16,15,13,14}, };

        //    // orismos shmeiwn
        //    for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
        //    {
        //        model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

        //    }

        //    // orismos elements 
        //    Element e1;
        //    int subdomainID = 1;
        //    for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
        //    {
        //        e1 = new Element()
        //        {
        //            ID = nElement + 1,
        //            ElementType = new Hexa8NLRAM_1(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
        //        };
        //        for (int j = 0; j < 8; j++)
        //        {
        //            e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
        //        }
        //        model.ElementsDictionary.Add(e1.ID, e1);
        //        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    }

        //    // constraint vashh opou z=-1
        //    for (int k = 1; k < 5; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }

        //    // fortish korufhs
        //    Load load1;
        //    for (int k = 17; k < 21; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.X,
        //            Amount = 1 * load_value
        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        //public static void HexaCantileverBuilder_copyMS(Model model, double load_value)
        //{
        //    ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };

        //    double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
        //    {0.250000,-0.250000,-1.000000},
        //    {-0.250000,0.250000,-1.000000},
        //    {0.250000,0.250000,-1.000000},
        //    {-0.250000,-0.250000,-0.500000},
        //    {0.250000,-0.250000,-0.500000},
        //    {-0.250000,0.250000,-0.500000},
        //    {0.250000,0.250000,-0.500000},
        //    {-0.250000,-0.250000,0.000000},
        //    {0.250000,-0.250000,0.000000},
        //    {-0.250000,0.250000,0.000000},
        //    {0.250000,0.250000,0.000000},
        //    {-0.250000,-0.250000,0.500000},
        //    {0.250000,-0.250000,0.500000},
        //    {-0.250000,0.250000,0.500000},
        //    {0.250000,0.250000,0.500000},
        //    {-0.250000,-0.250000,1.000000},
        //    {0.250000,-0.250000,1.000000},
        //    {-0.250000,0.250000,1.000000},
        //    {0.250000,0.250000,1.000000}};

        //    int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
        //    {2,12,11,9,10,8,7,5,6},
        //    {3,16,15,13,14,12,11,9,10},
        //    {4,20,19,17,18,16,15,13,14}, };

        //    // orismos shmeiwn
        //    for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
        //    {
        //        model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

        //    }

        //    // orismos elements 
        //    Element e1;
        //    int subdomainID = 1;
        //    for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
        //    {
        //        e1 = new Element()
        //        {
        //            ID = nElement + 1,
        //            ElementType = new Hexa8NLRAM_1copyMS(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
        //        };
        //        for (int j = 0; j < 8; j++)
        //        {
        //            e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
        //        }
        //        model.ElementsDictionary.Add(e1.ID, e1);
        //        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    }

        //    // constraint vashh opou z=-1
        //    for (int k = 1; k < 5; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }

        //    // fortish korufhs
        //    Load load1;
        //    for (int k = 17; k < 21; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.X,
        //            Amount = 1 * load_value
        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        #endregion

        public static void Hexa_1mat_CantileverBuilder(Model model, double load_value)
        {
            //xrhsimopoiithike to  ParadeigmataElegxwnBuilder.HexaCantileverBuilder(Model model, double load_value)
            // allagh tou element kai tou material

            //ElasticMaterial3D material1 = new ElasticMaterial3D()
            //{
            //    YoungModulus = 1353000,
            //    PoissonRatio = 0.3,
            //};


            VonMisesMaterial3D material1 = new VonMisesMaterial3D(1353000, 0.30, 1353000, 0.15);
            //{
            //    youngModulus = 1353000,
            //    PoissonRatio = 0.30,
            //};

            double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
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

            int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
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
            int subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NLRAM_1mat(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8                    
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
                model.NodesDictionary[k].Constraints.Add(DOFType.X);
                model.NodesDictionary[k].Constraints.Add(DOFType.Y);
                model.NodesDictionary[k].Constraints.Add(DOFType.Z);
            }

            // fortish korufhs
            Load load1;
            for (int k = 17; k < 21; k++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[k],
                    DOF = DOFType.X,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }

        public static void HexaElementsOnlyVonMises(Model model)
        {
            //proelefsi EmbeddedexamplesBuilder ekdoshs msolve telikhs pro feat/prosthiki_allagwn

            int startX = 0;
            int startY = 0;
            int startZ = 0;

            int nodeID = 1;
            for (int l = 0; l < 3; l++)
            {
                for (int k = 0; k < 2; k++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + j * 1, Y = startY + k * 1, Z = startZ + l * 1 });

                        nodeID++;
                    }
                }
            }
            nodeID = 1;
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    nodeID++;
                }
            }
            VonMisesMaterial3D material1 = new VonMisesMaterial3D(2.1e5, 0.35, 2.1e5, 0.15);
            //{
            //    youngModulus = 2.1e5,
            //    PoissonRatio = 0.35,
            //};

            //eisagwgh enos element
            Element e1;
            int ID2 = 1;
            e1 = new Element()
            {
                ID = 1,
                ElementType = new Hexa8(material1) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
            };
            ID2 = 1;
            for (int j = 0; j < 2; j++)
            {
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 1, model.NodesDictionary[4 * (ID2 - 1) + 1]); // na allaxthei h arithmisi swsth seira
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 2, model.NodesDictionary[4 * (ID2 - 1) + 2]);
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 4, model.NodesDictionary[4 * (ID2 - 1) + 4]);
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 3, model.NodesDictionary[4 * (ID2 - 1) + 3]);
                ID2++;
            }


            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            // ews edw

            //eisagwgh defterou element
            Element e2 = new Element()
            {
                ID = 2,
                ElementType = new Hexa8(material1) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
            };
            ID2 = 1;
            for (int j = 0; j < 2; j++)
            {
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 5, model.NodesDictionary[4 * (ID2 - 1) + 5]); // na allaxthei h arithmisi swsth seira
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 6, model.NodesDictionary[4 * (ID2 - 1) + 6]);
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 8, model.NodesDictionary[4 * (ID2 - 1) + 8]);
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 7, model.NodesDictionary[4 * (ID2 - 1) + 7]);
                ID2++;
            }
            model.ElementsDictionary.Add(e2.ID, e2);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
            // ews edw




            //model.Loads.Add()
            ID2 = 9;
            //apait 1
            DOFType doftype1;
            doftype1 = new DOFType();


            // apait 2
            Load load1;

            for (int j = 0; j < 4; j++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[ID2],
                    //DOF = doftype1,
                    DOF = DOFType.Z,
                    Amount = 500 //3*500// 2*500 //3*500 //250 //500

                };
                model.Loads.Add(load1);
                ID2++;
            }
        }

        //public static void HexaCantileverBuilderRAM(Model model, double load_value)
        //{
        //    ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };

        //    double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
        //    {0.250000,-0.250000,-1.000000},
        //    {-0.250000,0.250000,-1.000000},
        //    {0.250000,0.250000,-1.000000},
        //    {-0.250000,-0.250000,-0.500000},
        //    {0.250000,-0.250000,-0.500000},
        //    {-0.250000,0.250000,-0.500000},
        //    {0.250000,0.250000,-0.500000},
        //    {-0.250000,-0.250000,0.000000},
        //    {0.250000,-0.250000,0.000000},
        //    {-0.250000,0.250000,0.000000},
        //    {0.250000,0.250000,0.000000},
        //    {-0.250000,-0.250000,0.500000},
        //    {0.250000,-0.250000,0.500000},
        //    {-0.250000,0.250000,0.500000},
        //    {0.250000,0.250000,0.500000},
        //    {-0.250000,-0.250000,1.000000},
        //    {0.250000,-0.250000,1.000000},
        //    {-0.250000,0.250000,1.000000},
        //    {0.250000,0.250000,1.000000}};

        //    int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
        //    {2,12,11,9,10,8,7,5,6},
        //    {3,16,15,13,14,12,11,9,10},
        //    {4,20,19,17,18,16,15,13,14}, };

        //    // orismos shmeiwn
        //    for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
        //    {
        //        model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

        //    }

        //    // orismos elements 
        //    Element e1;
        //    int subdomainID = 1;
        //    for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
        //    {
        //        e1 = new Element()
        //        {
        //            ID = nElement + 1,
        //            ElementType = new Hexa8NLRAM(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
        //        };
        //        for (int j = 0; j < 8; j++)
        //        {
        //            e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
        //        }
        //        model.ElementsDictionary.Add(e1.ID, e1);
        //        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    }

        //    // constraint vashh opou z=-1
        //    for (int k = 1; k < 5; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }

        //    // fortish korufhs
        //    Load load1;
        //    for (int k = 17; k < 21; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.X,
        //            Amount = 1 * load_value
        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        //public static void HexaCantileverBuilderRAM2(Model model, double load_value)
        //{
        //    ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };

        //    double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
        //    {0.250000,-0.250000,-1.000000},
        //    {-0.250000,0.250000,-1.000000},
        //    {0.250000,0.250000,-1.000000},
        //    {-0.250000,-0.250000,-0.500000},
        //    {0.250000,-0.250000,-0.500000},
        //    {-0.250000,0.250000,-0.500000},
        //    {0.250000,0.250000,-0.500000},
        //    {-0.250000,-0.250000,0.000000},
        //    {0.250000,-0.250000,0.000000},
        //    {-0.250000,0.250000,0.000000},
        //    {0.250000,0.250000,0.000000},
        //    {-0.250000,-0.250000,0.500000},
        //    {0.250000,-0.250000,0.500000},
        //    {-0.250000,0.250000,0.500000},
        //    {0.250000,0.250000,0.500000},
        //    {-0.250000,-0.250000,1.000000},
        //    {0.250000,-0.250000,1.000000},
        //    {-0.250000,0.250000,1.000000},
        //    {0.250000,0.250000,1.000000}};

        //    int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
        //    {2,12,11,9,10,8,7,5,6},
        //    {3,16,15,13,14,12,11,9,10},
        //    {4,20,19,17,18,16,15,13,14}, };

        //    // orismos shmeiwn
        //    for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
        //    {
        //        model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

        //    }

        //    // orismos elements 
        //    Element e1;
        //    int subdomainID = 1;
        //    for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
        //    {
        //        e1 = new Element()
        //        {
        //            ID = nElement + 1,
        //            ElementType = new Hexa8NLRAM2(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
        //        };
        //        for (int j = 0; j < 8; j++)
        //        {
        //            e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
        //        }
        //        model.ElementsDictionary.Add(e1.ID, e1);
        //        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    }

        //    // constraint vashh opou z=-1
        //    for (int k = 1; k < 5; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }

        //    // fortish korufhs
        //    Load load1;
        //    for (int k = 17; k < 21; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.X,
        //            Amount = 1 * load_value
        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        //public static void HexaCantileverBuilderRAM_1(Model model, double load_value)
        //{
        //    ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };

        //    double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
        //    {0.250000,-0.250000,-1.000000},
        //    {-0.250000,0.250000,-1.000000},
        //    {0.250000,0.250000,-1.000000},
        //    {-0.250000,-0.250000,-0.500000},
        //    {0.250000,-0.250000,-0.500000},
        //    {-0.250000,0.250000,-0.500000},
        //    {0.250000,0.250000,-0.500000},
        //    {-0.250000,-0.250000,0.000000},
        //    {0.250000,-0.250000,0.000000},
        //    {-0.250000,0.250000,0.000000},
        //    {0.250000,0.250000,0.000000},
        //    {-0.250000,-0.250000,0.500000},
        //    {0.250000,-0.250000,0.500000},
        //    {-0.250000,0.250000,0.500000},
        //    {0.250000,0.250000,0.500000},
        //    {-0.250000,-0.250000,1.000000},
        //    {0.250000,-0.250000,1.000000},
        //    {-0.250000,0.250000,1.000000},
        //    {0.250000,0.250000,1.000000}};

        //    int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
        //    {2,12,11,9,10,8,7,5,6},
        //    {3,16,15,13,14,12,11,9,10},
        //    {4,20,19,17,18,16,15,13,14}, };

        //    // orismos shmeiwn
        //    for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
        //    {
        //        model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

        //    }

        //    // orismos elements 
        //    Element e1;
        //    int subdomainID = 1;
        //    for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
        //    {
        //        e1 = new Element()
        //        {
        //            ID = nElement + 1,
        //            ElementType = new Hexa8NLRAM_1(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
        //        };
        //        for (int j = 0; j < 8; j++)
        //        {
        //            e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
        //        }
        //        model.ElementsDictionary.Add(e1.ID, e1);
        //        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    }

        //    // constraint vashh opou z=-1
        //    for (int k = 1; k < 5; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }

        //    // fortish korufhs
        //    Load load1;
        //    for (int k = 17; k < 21; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.X,
        //            Amount = 1 * load_value
        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        //public static void Example2Hexa8NL1Cohesive8node(Model model)
        //{
        //    ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 135300, // 1353000 gia to allo paradeigma
        //        PoissonRatio = 0.3,
        //    };
        //    BenzeggaghKenaneCohMat material2 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57,// N / mm2
        //        D_o_3 = 0.000057, // mm
        //        D_f_3 = 0.0098245610,  // mm

        //        T_o_1 = 57,// N / mm2
        //        D_o_1 = 0.000057, // mm
        //        D_f_1 = 0.0098245610,  // mm

        //        n_curve = 1.4
        //    };


        //    double[,] nodeData = new double[,] {
        //    {0.500000,0.000000,1.000000},
        //    {0.500000,0.500000,1.000000},
        //    {0.000000,0.000000,1.000000},
        //    {0.000000,0.500000,1.000000},
        //    {0.500000,0.000000,0.500000},
        //    {0.500000,0.500000,0.500000},
        //    {0.000000,0.000000,0.500000},
        //    {0.000000,0.500000,0.500000},
        //    {0.500000,0.000000,0.000000},
        //    {0.500000,0.500000,0.000000},
        //    {0.000000,0.000000,0.000000},
        //    {0.000000,0.500000,0.000000},
        //    {0.500000,0.000000,0.500000},
        //    {0.500000,0.500000,0.500000},
        //    {0.000000,0.000000,0.500000},
        //    {0.000000,0.500000,0.500000} };

        //    int[,] elementData = new int[,] {{1,1,2,4,3,5,6,8,7},{2,13,14,16,15,9,10,12,11},};

        //    int[] cohesive8Nodes = new int[] { 5, 6, 8, 7, 13, 14, 16, 15, };

        //    // orismos shmeiwn
        //    for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
        //    {
        //        model.NodesDictionary.Add(nNode+1, new Node() { ID = nNode+1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

        //    }

        //    // orismos elements 
        //    Element e1;
        //    int subdomainID = 1;
        //    for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
        //    {
        //        e1 = new Element()
        //        {
        //            ID = nElement+1,
        //            ElementType = new Hexa8NL(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
        //        };
        //        for (int j = 0; j < 8; j++)
        //        {
        //            e1.NodesDictionary.Add(elementData[nElement,j+1], model.NodesDictionary[elementData[nElement, j+1]]);
        //        }
        //        model.ElementsDictionary.Add(e1.ID, e1);
        //        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    }

        //    // kai to cohesive
        //    e1 = new Element()
        //    {
        //        ID = 3,
        //        ElementType = new cohesive8node(material2, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
        //    };
        //    for (int j = 0; j < 8; j++)
        //    {
        //        e1.NodesDictionary.Add(cohesive8Nodes[j], model.NodesDictionary[cohesive8Nodes[j]]);
        //    }
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //}

        //public static void Example2Hexa8NL1Cohesive8nodeConstraintsOrthiElastic(Model model)
        //{
        //    for (int k = 13; k < 17; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }
        //}

        //public static void Example2Hexa8NL1Cohesive8nodeConstraintsMixed(Model model)
        //{            
        //    for (int k = 13; k < 17; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }
        //}

        //public static void Example2Hexa8NL1Cohesive8nodeConstraintsBenc1(Model model)
        //{
        //    for (int k = 1; k < 5; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }
        //    for (int k = 5; k < 9; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //    }
        //    for (int k = 9; k < 17; k++)
        //    {
        //        model.NodesDictionary[k].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[k].Constraints.Add(DOFType.Z);
        //    }
        //}

        //public static void Example2Hexa8NL1Cohesive8nodeLoadsMIxed(Model model, double load_value)
        //{
        //    Load load1;
        //    for (int k = 5; k < 9; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.X,
        //            Amount = 1*load_value

        //        };
        //        model.Loads.Add(load1);
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.Z,
        //            Amount = 0.2*load_value

        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        //public static void Example2Hexa8NL1Cohesive8nodeLoadsElasticOrthi(Model model, double load_value)
        //{
        //    Load load1;
        //    for (int k = 1; k < 5; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.Z,
        //            Amount = 1 * load_value

        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        //public static void Example2Hexa8NL1Cohesive8nodeLoadsBenc1(Model model, double load_value)
        //{
        //    Load load1;
        //    for (int k = 5; k < 9; k++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[k],
        //            DOF = DOFType.Z,
        //            Amount =1* load_value

        //        };
        //        model.Loads.Add(load1);
        //    }
        //}

        # region tocomment out temporarily 
        //public static void ShellAndCohesiveShellPaktwsh(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = 1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57, // New load case argurhs NR_shell_coh.m
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D_v2 material2 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlk(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2 = 0,
        //        }
        //    };

        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2*2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}
        #endregion

        //public static void ShellAndCohesiveRAM_1ShellPaktwsh(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = 1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57, // New load case argurhs NR_shell_coh.m
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D_v2 material2 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2 = 0,
        //        }
        //    };

        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2 * 2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}

        //public static void ShellAndCohesiveRAM_11ShellPaktwsh(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = 1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57, // New load case argurhs NR_shell_coh.m
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D_v2 material2 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_11(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2 = 0,
        //        }
        //    };

        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2 * 2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}

        //public static void ShellAndCohesiveShellPaktwshRAM(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = 1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57, // New load case argurhs NR_shell_coh.m
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D_v2 material2 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyGetRAM(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGet(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2 = 0,
        //        }
        //    };

        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2 * 2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}

        //public static void ShellAndCohesiveShellPaktwshRAM2(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = 1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57, // New load case argurhs NR_shell_coh.m
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D_v2 material2 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyGetRAM2(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGet(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2 = 0,
        //        }
        //    };

        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2 * 2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}

        //public static void ShellAndCohesiveRAMShellPaktwsh(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = 1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57, // New load case argurhs NR_shell_coh.m
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D_v2 material2 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2 = 0,
        //        }
        //    };

        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2 * 2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}

        //public static void ShellAndCohesiveShellRAM_1Paktwsh(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = 1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57, // New load case argurhs NR_shell_coh.m
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D_v2 material2 = new ElasticMaterial3D_v2()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2 = 0,
        //        }
        //    };

        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2 * 2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}

    }
}
