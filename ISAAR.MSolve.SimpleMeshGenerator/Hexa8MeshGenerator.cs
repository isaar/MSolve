using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.SimpleMeshGenerator
{
    public enum Surface
    {
        None,
	    XZDown, 
        XZUp,
	    YZLeft, 
        YZRight,
	    XYFront, 
        XYBack,
        All
    }

    public static class Hexa8MeshGenerator
    {
        private static int[] xyzOrder = new int[] { 0, 1, 2 };
        private static int[] nodesPerElementSide = new int[] { 2, 2, 2 };
        private static readonly Hexa8Memoizer hexa8Memoizer = new Hexa8Memoizer();
        public static readonly Model model = new Model();

        private static void PinSurface(int[] nodesPerDirection, Surface pinSurface, int nodeID, int i, int j, int k, bool hasDampers)
        {
            if (hasDampers)
            {
                if (k == 0)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                }
                if (k == nodesPerDirection[0] + 1)
                //if (i == nodesPerDirection[0] - 1)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                }
        }

            switch (pinSurface)
            {
                case Surface.XZDown:
                    if (j == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    break;
                case Surface.XZUp:
                    if (j == nodesPerDirection[1] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    break;
                case Surface.YZLeft:
                    if (i == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    break;
                case Surface.YZRight:
                    if (k == nodesPerDirection[0] - 1)
                    //if (i == nodesPerDirection[0] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    break;
                case Surface.XYFront:
                    if (k == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    break;
                case Surface.XYBack:
                    if (k == nodesPerDirection[2] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    break;
                case Surface.All:
                    if (j == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    if (j == nodesPerDirection[1] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    if (i == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    if (i == nodesPerDirection[0] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    if (k == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    if (k == nodesPerDirection[2] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                    break;
            }
        }

        private static void PoreSurface(int[] nodesPerDirection, Surface poreSurface, int nodeID, int i, int j, int k)
        {
            switch (poreSurface)
            {
                case Surface.XZDown:
                    if (j == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Pore);
                    }
                    break;
                case Surface.XZUp:
                    if (j == nodesPerDirection[1] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Pore);
                    }
                    break;
                case Surface.YZLeft:
                    if (i == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Pore);
                    }
                    break;
                case Surface.YZRight:
                    if (i == nodesPerDirection[0] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Pore);
                    }
                    break;
                case Surface.XYFront:
                    if (k == 0)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Pore);
                    }
                    break;
                case Surface.XYBack:
                    if (k == nodesPerDirection[2] - 1)
                    {
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Pore);
                    }
                    break;
            }
        }

        private static void LoadSurface(int[] nodesPerDirection, Surface loadSurface, int nodeID, 
            DOFType loadDirection, double loadAmount, int i, int j, int k)
        {
            switch (loadSurface)
            {
                case Surface.None:
                    if (j == nodesPerDirection[1] - 1 && k == nodesPerDirection[0] - 1)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
                case Surface.XZDown:
                    if (j == 0)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
                case Surface.XZUp:
                    if (j == nodesPerDirection[1] - 1)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
                case Surface.YZLeft:
                    if (i == 0)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
                case Surface.YZRight:
                    if (i == nodesPerDirection[0] - 1)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
                case Surface.XYFront:
                    if (k == 0)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
                case Surface.XYBack:
                    if (k == nodesPerDirection[2] - 1)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
                case Surface.All:
                    if (j == 0)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    if (j == nodesPerDirection[1] - 1)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    if (i == 0)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    if (i == nodesPerDirection[0] - 1)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    if (k == 0)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    if (k == nodesPerDirection[2] - 1)
                    {
                        model.Loads.Add(new Load()
                        {
                            Node = model.NodesDictionary[nodeID],
                            DOF = loadDirection,
                            Amount = loadAmount
                        });
                    }
                    break;
            }
        }

        private static void LoadSurfaceMassAcceleration(int[] elementsPerDirection, Surface loadSurface, int elementID,
            MassAccelerationHistoryLoad historyLoad, int i, int j, int k)
        {
            switch (loadSurface)
            {
                case Surface.XZDown:
                    if (j == 0)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    break;
                case Surface.XZUp:
                    if (j == elementsPerDirection[1] - 1)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    break;
                case Surface.YZLeft:
                    if (i == 0)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    break;
                case Surface.YZRight:
                    if (i == elementsPerDirection[0] - 1)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    break;
                case Surface.XYFront:
                    if (k == 0)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    break;
                case Surface.XYBack:
                    if (k == elementsPerDirection[2] - 1)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    break;
                case Surface.All:
                    if (j == 0)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    if (j == elementsPerDirection[1] - 1)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    if (i == 0)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    if (i == elementsPerDirection[0] - 1)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    if (k == 0)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    if (k == elementsPerDirection[2] - 1)
                    {
                        model.ElementMassAccelerationHistoryLoads.Add(new ElementMassAccelerationHistoryLoad()
                        {
                            Element = model.ElementsDictionary[elementID],
                            HistoryLoad = historyLoad
                        });
                    }
                    break;
            }
        }

        private static void BuildConstraintsAndLoads(int[] nodesPerDirection, Surface pinSurface, Surface poreSurface, 
            Surface loadSurface, DOFType loadDirection, double loadAmount, bool hasDampers)
        {
            int nodeID = 0;
            if (hasDampers)
            {
                for (int j = 0; j < nodesPerDirection[xyzOrder[1]]; j++)
                {
                    for (int k = 0; k < nodesPerDirection[xyzOrder[0]]; k++)
                    {
                        nodeID++;
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z); 
                    }
                }
            }

            for (int i = 0; i < nodesPerDirection[2]; i++)
            {
                for (int j = 0; j < nodesPerDirection[1]; j++)
                {
                    for (int k = 0; k < nodesPerDirection[0] + (hasDampers ? 2 : 0); k++)
                    {
                        nodeID++;
                        PoreSurface(nodesPerDirection, poreSurface, nodeID, i, j, k);
                        PinSurface(nodesPerDirection, pinSurface, nodeID, i, j, k, hasDampers);
                        LoadSurface(nodesPerDirection, loadSurface, nodeID, loadDirection, loadAmount, i, j, k);
                    }
                }
            }

            if (hasDampers)
            {
                for (int j = 0; j < nodesPerDirection[xyzOrder[1]]; j++)
                {
                    for (int k = 0; k < nodesPerDirection[xyzOrder[0]]; k++)
                    {
                        nodeID++;
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    }
                }
            }
        }

        private static void BuildElementMassAccelerationLoads(int[] elementsPerDirection, Surface loadSurface, 
            MassAccelerationHistoryLoad loadAmount)
        {
            int elementID = 0;
            for (int i = 0; i < elementsPerDirection[2]; i++)
            {
                for (int j = 0; j < elementsPerDirection[1]; j++)
                {
                    for (int k = 0; k < elementsPerDirection[0]; k++)
                    {
                        elementID++;
                        LoadSurfaceMassAcceleration(elementsPerDirection, loadSurface, elementID, loadAmount, i, j, k);
                    }
                }
            }
        }

        private static void BuildNodes(int[] nodesPerDirection, double[] distancesPerDirection, bool hasDampers)
        {
            int nodes = 0;
            double[] xyz = new double[] { hasDampers ? distancesPerDirection[xyzOrder[0]] : 0, 0, 0 };
            if (hasDampers)
            {
                for (int j = 0; j < nodesPerDirection[xyzOrder[1]]; j++)
                {
                    for (int k = 0; k < nodesPerDirection[xyzOrder[0]]; k++)
                    {
                        nodes++;
                        Node node = new Node() { ID = nodes, X = xyz[0], Y = xyz[1], Z = xyz[2] };
                        model.NodesDictionary.Add(nodes, node);
                        xyz[xyzOrder[0]] += distancesPerDirection[xyzOrder[0]];
                    }
                    xyz[xyzOrder[0]] = distancesPerDirection[xyzOrder[0]];
                    xyz[xyzOrder[1]] += distancesPerDirection[xyzOrder[1]];
                }
                xyz[xyzOrder[0]] = 0;
                xyz[xyzOrder[1]] = 0;
                xyz[xyzOrder[2]] += distancesPerDirection[xyzOrder[2]];
            }

            for (int i = 0; i < nodesPerDirection[xyzOrder[2]]; i++)
            {
                for (int j = 0; j < nodesPerDirection[xyzOrder[1]]; j++)
                {
                    for (int k = 0; k < nodesPerDirection[xyzOrder[0]] + (hasDampers ? 2 : 0); k++)
                    {
                        nodes++;
                        Node node = new Node() { ID = nodes, X = xyz[0], Y = xyz[1], Z = xyz[2] };
                        model.NodesDictionary.Add(nodes, node);
                        xyz[xyzOrder[0]] += distancesPerDirection[xyzOrder[0]];
                    }
                    xyz[xyzOrder[0]] = 0;
                    xyz[xyzOrder[1]] += distancesPerDirection[xyzOrder[1]];
                }
                xyz[xyzOrder[1]] = 0;
                xyz[xyzOrder[2]] += distancesPerDirection[xyzOrder[2]];
            }

            if (hasDampers)
            {
                xyz[xyzOrder[0]] = distancesPerDirection[xyzOrder[0]];
                for (int j = 0; j < nodesPerDirection[xyzOrder[1]]; j++)
                {
                    for (int k = 0; k < nodesPerDirection[xyzOrder[0]]; k++)
                    {
                        nodes++;
                        Node node = new Node() { ID = nodes, X = xyz[0], Y = xyz[1], Z = xyz[2] };
                        model.NodesDictionary.Add(nodes, node);
                        xyz[xyzOrder[0]] += distancesPerDirection[xyzOrder[0]];
                    }
                    xyz[xyzOrder[0]] = distancesPerDirection[xyzOrder[0]];
                    xyz[xyzOrder[1]] += distancesPerDirection[xyzOrder[1]];
                }
                xyz[xyzOrder[0]] = 0;
                xyz[xyzOrder[1]] = 0;
                xyz[xyzOrder[2]] += distancesPerDirection[xyzOrder[2]];
            }
        }

        //private static IFiniteElementMaterial3D GetElasticMaterial()
        //{
        //    return new ElasticMaterial3D()
        //    {
        //        YoungModulus = 250000,
        //        PoissonRatio = 0.3
        //    };
        //}

        private static IStochasticFiniteElementMaterial GetStochasticPlasticMaterial(IStochasticMaterialCoefficientsProvider p, int pos = 0)
        {
            return new StochasticElasticMaterial3D(p);
            //switch (pos)
            //{
            //    case 1:
            //        return new StochasticMohrCoulombMaterial(p, 6000, 0.25, 17.5, 17, 0);
            //    case 2:
            //        return new StochasticMohrCoulombMaterial(p, 20000, 0.3, 0, 35, 5);
            //    case 3:
            //        return new StochasticMohrCoulombMaterial(p, 60000, 0.35, 300, 40, 0);
            //    default:
            //        return new StochasticMohrCoulombMaterial(p, 60000, 0.35, 300, 40, 0);
            //}
        }

        private static IStochasticFiniteElementMaterial GetStochasticElasticMaterial(IStochasticMaterialCoefficientsProvider p)
        {
            return new StochasticElasticMaterial3D(p)
            {
                YoungModulus = 21000000, // New Giovanis
                //YoungModulus = 125000000, // Aris
                //YoungModulus = 1000000, (Giovanis)
                PoissonRatio = 0.3
            };
        }

        private static IStochasticFiniteElementMaterial GetStifferStochasticElasticMaterial(IStochasticMaterialCoefficientsProvider p)
        {
            return new StochasticElasticMaterial3D(p)
            {
                YoungModulus = 2500000000,
                PoissonRatio = 0.3
            };
        }

        private static IFiniteElementMaterial3D GetHostElasticMaterial()
        {
            return new ElasticMaterial3D()
            {
                //YoungModulus = 3.1E+10,
                //PoissonRatio = 0.2
                YoungModulus = 2.8,
                PoissonRatio = 0.4
            };
        }

        private static IFiniteElementMaterial3D GetEmbeddedElasticMaterial()
        {
            return new ElasticMaterial3D()
            {
                //YoungModulus = 2.1E+11,
                //PoissonRatio = 0.32
                YoungModulus = 1051,
                PoissonRatio = 0.04473161
            };
        }

        private static IFiniteElementMaterial3D GetElasticMaterial(int pos = 0)
        {
            double E = 0;
            double n = 0;
            switch (pos)
            {
                case 1:
                    E = 6000;
                    n = 0.25;
                    break;
                case 2:
                    E = 20000;
                    n = 0.3;
                    break;
                case 3:
                    E = 60000;
                    n = 0.35;
                    break;
                default:
                    //E = 100000, // Ambrosios
                    //E = 21000000, // New Giovanis
                    E = 100000; // Goat tet
                    //E = 125000000, // Aris
                    //E = 1000000, (Giovanis)
                    n = 0.3;
                    break;
            }
            return new ElasticMaterial3D()
            {
                YoungModulus = E, 
                PoissonRatio = n
            };
        }

        private static IFiniteElementMaterial3D GetStifferElasticMaterial()
        {
            return new ElasticMaterial3D()
            {
                YoungModulus = 2500000000,
                PoissonRatio = 0.3
            };
        }

        private static IFiniteElementMaterial3D GetPlasticMaterial(int pos = 0)
        {
            //switch (pos)
            //{
            //    case 1:
            //        return new MohrCoulombMaterial(12000, 0.25, 16.2, 30, 0);
            //        //return new VonMisesMaterial3D(12000, 0.25, 52, 0);
            //    case 2:
            //        return new MohrCoulombMaterial(50000, 0.3, 90, 20, 0);
            //        //return new VonMisesMaterial3D(50000, 0.3, 125, 0);
            //    case 3:
            //        return new MohrCoulombMaterial(150000, 0.35, 250, 38, 0);
            //        //return new VonMisesMaterial3D(150000, 0.35, 350, 0);
            //    default:
            //        return new VonMisesMaterial3D(100000, 0.3, 100, 0);
            //}
            switch (pos)
            {
                case 1:
                    return new MohrCoulombMaterial(6000, 0.25, 17.5, 17, 0);
                //return new VonMisesMaterial3D(12000, 0.25, 52, 0);
                case 2:
                    return new MohrCoulombMaterial(20000, 0.3, 0, 35, 5);
                //return new VonMisesMaterial3D(50000, 0.3, 125, 0);
                case 3:
                    return new MohrCoulombMaterial(60000, 0.35, 300, 40, 0);
                //return new VonMisesMaterial3D(150000, 0.35, 350, 0);
                default:
                    return new VonMisesMaterial3D(1000000, 0.3, 10, 0);
            }

            //return new VonMisesMaterial3DFORTRAN()
            //{
            //    YoungModulus = 250000,
            //    PoissonRatio = 0.3,
            //    YieldStress = 255,
            //    HardeningRatio = 3900
            //};
        }

        private static IFiniteElement GetSolidElementType(bool isElastic, bool isStiff, bool isStochastic, bool isHost, IStochasticMaterialCoefficientsProvider p, int pos = 0)
        {
            double porosity = 0;
            double density = 0;
            switch (pos)
            {
                case 1:
                    porosity = 0.350;
                    density = 1.9;
                    break;
                case 2:
                    porosity = 0.455;
                    density = 1.8;
                    break;
                case 3:
                    porosity = 0.3;
                    density = 1.85;
                    break;
                default:
                    porosity = 0.455;
                    density = 2.2;
                    break;
            }

            if (isHost) return new Hexa8(GetHostElasticMaterial()) { Density = porosity * 1 * 1 + (1 - porosity) * density };

            //if (p is LognormalPCFileStochasticCoefficientsProvider || p is GaussianPCFileStochasticCoefficientsProvider)
            //    return isStiff ?
            //        new Hexa8Stochastic(GetStifferElasticMaterial()) { Density = porosity * 1 * 1 + (1 - porosity) * density } :
            //        new Hexa8Stochastic(GetElasticMaterial(pos)) { Density = porosity * 1 * 1 + (1 - porosity) * density };

            if (p != null)
                if (isElastic)
                    return isStiff ?
                        new Hexa8WithStochasticMaterial(GetStifferStochasticElasticMaterial(p), hexa8Memoizer) { Density = porosity * 1 * 1 + (1 - porosity) * density } :
                        new Hexa8WithStochasticMaterial(GetStochasticElasticMaterial(p), hexa8Memoizer) { Density = porosity * 1 * 1 + (1 - porosity) * density };
                else
                    return
                        new Hexa8WithStochasticMaterial(GetStochasticPlasticMaterial(p, pos), hexa8Memoizer) { Density = porosity * 1 * 1 + (1 - porosity) * density };

            IFiniteElementMaterial3D m;
            if (isElastic)
                if (isStiff)
                    m = GetStifferElasticMaterial();
                else
                    m = GetElasticMaterial(pos);
            else
                m = GetPlasticMaterial(pos);

            ////return isStochastic ? new Hexa8Stochastic(m) { Density = 0.455 * 1 * 1 + (1 - 0.455) * 2.2 } :
            //return isStochastic ? new Hexa8Stochastic(m) { Density = porosity * 1 * 1 + (1 - porosity) * density } :
            //    new Hexa8(m)
            //    {
            //        Density = porosity * 1 * 1 + (1 - porosity) * density,
            //        //RayleighAlpha = 0.347114713,
            //        //RayleighBeta = 1.5315E-05,
            //    };
            ////            new Hexa8(m) { Density = 0.455 * 1 * 1 + (1 - 0.455) * 2.2 };

            return new Hexa8(m)
            {
                Density = porosity * 1 * 1 + (1 - porosity) * density,
                //RayleighAlpha = 0.347114713,
                //RayleighBeta = 1.5315E-05,
            };
        }

        private static IFiniteElement GetPorousElementType(bool isElastic, IStochasticMaterialCoefficientsProvider p, int pos = 0)
        {
            IFiniteElementMaterial3D m = null;
            IStochasticFiniteElementMaterial mm = null;
            if (p == null)
            {
                if (isElastic)
                    m = GetElasticMaterial(pos);
                else
                    m = GetPlasticMaterial(pos);
            }
            else
            {
                if (isElastic)
                    throw new ArgumentException("Stochastic porous with elastic material not supported.");
                else
                    mm = GetStochasticPlasticMaterial(p, pos);
            }
            double porosity = 0;
            double density = 0;
            double permeability = 0;
            //switch (pos)
            //{
            //    case 1:
            //        porosity = 0.35;
            //        density = 1.8;
            //        permeability = 1.019e-7;
            //        break;
            //    case 2:
            //        porosity = 0.45;
            //        density = 1.95;
            //        permeability = 3.058e-10;
            //        break;
            //    case 3:
            //        porosity = 0.25;
            //        density = 2.1;
            //        permeability = 1.019e-4;
            //        break;
            //    default:
            //        porosity = 0.455;
            //        density = 2.2;
            //        permeability = 0.000000012;
            //        break;
            //}
            switch (pos)
            {
                case 1:
                    porosity = 0.35;
                    density = 1.9;
                    permeability = 1.019e-7;
                    break;
                case 2:
                    porosity = 0.455;
                    density = 1.8;
                    permeability = 3.058e-4;
                    break;
                case 3:
                    porosity = 0.3;
                    density = 1.85;
                    permeability = 1.019e-5;
                    break;
                default:
                    porosity = 0.455;
                    density = 2.2;
                    permeability = 0.000000012;
                    break;
            }

            return p == null ? new Hexa8u8p(m) 
            {
                //RayleighAlpha = 4.313748979,
                //RayleighBeta = 0.013612851,
                //RayleighAlpha = 0.345099918,
                //RayleighBeta = 0.001089028,
                //RayleighAlpha = 0.347114713,
                //RayleighBeta = 1.5315E-05,
                FluidBulkModulus = 2200000,
                FluidDensity = 1,
                Permeability = permeability,
                PoreA = 1,
                Porosity = porosity,
                Saturation = 1,
                Xw = 1,
                SolidDensity = porosity * 1 * 1 + (1 - porosity) * density
            } :
            new Hexa8u8pWithStochasticMaterial(mm)
            {
                FluidBulkModulus = 2200000,
                FluidDensity = 1,
                Permeability = permeability,
                PoreA = 1,
                Porosity = porosity,
                Saturation = 1,
                Xw = 1,
                SolidDensity = porosity * 1 * 1 + (1 - porosity) * density
            };
        }

        private static IFiniteElement GetSpringElementType(bool isX, int pos = 0)
        {
            var spring = 0d;
            var damper = 0d;
            //switch (pos)
            //{
            //    case 1:
            //        spring = 12000;
            //        damper = 2.75e3;
            //        break;
            //    case 2:
            //        spring = 50000;
            //        damper = 5.61e3;
            //        break;
            //    case 3:
            //        spring = 150000;
            //        damper = 9.72e3;
            //        break;
            //    default:
            //        spring = 100000;
            //        damper = 2.75e4;
            //        break;
            //}
            switch (pos)
            {
                case 1:
                    spring = 12000;
                    damper = 1.904256e3;
                    break;
                case 2:
                    spring = 40000;
                    damper = 3.621709e3;
                    break;
                case 3:
                    spring = 120000;
                    damper = 6.57763e3;
                    break;
                default:
                    spring = 100000;
                    damper = 2.75e4;
                    break;
            }
            return new SpringDamper3D(spring, damper, isX ? SpringDirections.X : SpringDirections.Z, SpringDirections.XYZ);
        }

        private static int[] BuildElements(int[] nodesPerDirection, bool isPorous, bool isElastic, bool isStochastic, bool isIllConditioned, IStochasticMaterialCoefficientsProvider p, bool hasDampers)
        {
            int[] elementsPerDirection = new int[] 
                { 
                    (int)((nodesPerDirection[0] - 1) / (nodesPerElementSide[0] - 1)),
                    (int)((nodesPerDirection[1] - 1) / (nodesPerElementSide[1] - 1)),
                    (int)((nodesPerDirection[2] - 1) / (nodesPerElementSide[2] - 1))
                };
            int nodeStart = 0;
            int elements = 0;
            int springElements = 0;
            int nX = nodesPerDirection[0] + (hasDampers ? 2 : 0);
            int nY = nodesPerDirection[1];
            int nZ = nodesPerDirection[2];

            if (hasDampers)
            {
                int nnX = nodesPerDirection[0];
                Element element;
                nodeStart = nodesPerDirection[xyzOrder[1]] * nodesPerDirection[xyzOrder[0]];
                for (int j = 0; j < nodesPerDirection[xyzOrder[1]]; j++)
                {
                    int elProp = 3;
                    if (j >= 1)
                        elProp = 2;
                    if (j >= 3)
                        elProp = 1;
                    for (int k = 0; k < nodesPerDirection[xyzOrder[0]]; k++)
                    {
                        springElements--;
                        element = new Element();
                        element.ElementType = GetSpringElementType(false, elProp);
                        element.ID = springElements;
                        element.AddNode(model.NodesDictionary[nnX * j + k + 1]);
                        element.AddNode(model.NodesDictionary[nodeStart + (j * 2) + 1 + nnX * j + k + 1]);
                        model.ElementsDictionary.Add(element.ID, element);
                    }
                }
            }

            for (int i = 0; i < elementsPerDirection[2]; i++)
                for (int j = 0; j < elementsPerDirection[1]; j++)
                    for (int k = 0; k < elementsPerDirection[0]; k++)
                    {
                        int elProp = 3;
                        if (j >= 1)
                            elProp = 2;
                        if (j >= 3)
                            elProp = 1;
                        
                        // Uncomment if not SSI example                        
                        //elProp = 0;
                        Element element;

                        elements++;
                        if (hasDampers)
                        {
                            if (k == 0)
                            {
                                springElements--;
                                element = new Element();
                                element.ElementType = GetSpringElementType(true, elProp);
                                element.ID = springElements;
                                element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * j + 1]);
                                element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * j + 2]);
                                model.ElementsDictionary.Add(element.ID, element);

                                if (i == elementsPerDirection[2] - 1)
                                {
                                    springElements--;
                                    element = new Element();
                                    element.ElementType = GetSpringElementType(true, elProp);
                                    element.ID = springElements;
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * j + 1]);
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * j + 2]);
                                    model.ElementsDictionary.Add(element.ID, element);
                                }

                                if (j == elementsPerDirection[1] - 1)
                                {
                                    springElements--;
                                    element = new Element();
                                    element.ElementType = GetSpringElementType(true, elProp);
                                    element.ID = springElements;
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * elementsPerDirection[1] + 1]);
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * elementsPerDirection[1] + 2]);
                                    model.ElementsDictionary.Add(element.ID, element);

                                    if (i == elementsPerDirection[2] - 1)
                                    {
                                        springElements--;
                                        element = new Element();
                                        element.ElementType = GetSpringElementType(true, elProp);
                                        element.ID = springElements;
                                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * elementsPerDirection[1] + 1]);
                                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * elementsPerDirection[1] + 2]);
                                        model.ElementsDictionary.Add(element.ID, element);
                                    }
                                }
                            }
                            else if (k == elementsPerDirection[0] - 1)
                            {
                                springElements--;
                                element = new Element();
                                element.ElementType = GetSpringElementType(true, elProp);
                                element.ID = springElements;
                                element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * j + nodesPerDirection[0] + 1]);
                                element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * j + nodesPerDirection[0] + 2]);
                                model.ElementsDictionary.Add(element.ID, element);

                                if (i == elementsPerDirection[2] - 1)
                                {
                                    springElements--;
                                    element = new Element();
                                    element.ElementType = GetSpringElementType(true, elProp);
                                    element.ID = springElements;
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * j + nodesPerDirection[0] + 1]);
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * j + nodesPerDirection[0] + 2]);
                                    model.ElementsDictionary.Add(element.ID, element);
                                }

                                if (j == elementsPerDirection[1] - 1)
                                {
                                    springElements--;
                                    element = new Element();
                                    element.ElementType = GetSpringElementType(true, elProp);
                                    element.ID = springElements;
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * elementsPerDirection[1] + nodesPerDirection[0] + 1]);
                                    element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * elementsPerDirection[1] + nodesPerDirection[0] + 2]);
                                    model.ElementsDictionary.Add(element.ID, element);

                                    if (i == elementsPerDirection[2] - 1)
                                    {
                                        springElements--;
                                        element = new Element();
                                        element.ElementType = GetSpringElementType(true, elProp);
                                        element.ID = springElements;
                                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * elementsPerDirection[1] + nodesPerDirection[0] + 1]);
                                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * elementsPerDirection[2] + nX * elementsPerDirection[1] + nodesPerDirection[0] + 2]);
                                        model.ElementsDictionary.Add(element.ID, element);
                                    }
                                }
                            }
                        }

                        element = new Element();
                        if (isPorous)
                            element.ElementType = GetPorousElementType(isElastic, p, elProp);
                        else
                            element.ElementType = GetSolidElementType(isElastic, false, isStochastic, false, p, elProp);
                        element.ID = elements;
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * j + k + 1 + (hasDampers ? 1 : 0)]);
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * j + k + 2 + (hasDampers ? 1 : 0)]);
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * (j + 1) + k + 2 + (hasDampers ? 1 : 0)]);
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * i + nX * (j + 1) + k + 1 + (hasDampers ? 1 : 0)]);
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * (i + 1) + nX * j + k + 1 + (hasDampers ? 1 : 0)]);
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * (i + 1) + nX * j + k + 2 + (hasDampers ? 1 : 0)]);
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * (i + 1) + nX * (j + 1) + k + 2 + (hasDampers ? 1 : 0)]);
                        element.AddNode(model.NodesDictionary[nodeStart + nX * nY * (i + 1) + nX * (j + 1) + k + 1 + (hasDampers ? 1 : 0)]);
                        model.ElementsDictionary.Add(element.ID, element);
                    }

            if (hasDampers)
            {
                int nodeStart2 = nodeStart + nX * nY * elementsPerDirection[2];
                int nnX = nodesPerDirection[0];
                Element element;
                nodeStart = nX * nY;
                for (int j = 0; j < nodesPerDirection[xyzOrder[1]]; j++)
                {
                    int elProp = 3;
                    if (j >= 1)
                        elProp = 2;
                    if (j >= 3)
                        elProp = 1;
                    for (int k = 0; k < nodesPerDirection[xyzOrder[0]]; k++)
                    {
                        springElements--;
                        element = new Element();
                        element.ElementType = GetSpringElementType(false, elProp);
                        element.ID = springElements;
                        element.AddNode(model.NodesDictionary[nodeStart2 + (j * 2) + 1 + nnX * j + k + 1]);
                        element.AddNode(model.NodesDictionary[nodeStart2 + nodeStart + nnX * j + k + 1]);
                        model.ElementsDictionary.Add(element.ID, element);
                    }
                }
            }
            if (isIllConditioned)
            {
                if (elementsPerDirection[0] != 10 || elementsPerDirection[1] != 10 || isPorous || !isElastic)
                    throw new ArgumentException("Ill conditioning is for a 10x10xN elastic solid cube.");

                int faceElements = elementsPerDirection[0] * elementsPerDirection[1];
                for (int i = 0; i < elementsPerDirection[2]; i++)
                {
                    model.ElementsDictionary[i * faceElements + 33].ElementType = GetSolidElementType(isElastic, true, isStochastic, false, p);
                    model.ElementsDictionary[i * faceElements + 38].ElementType = GetSolidElementType(isElastic, true, isStochastic, false, p);
                    model.ElementsDictionary[i * faceElements + 73].ElementType = GetSolidElementType(isElastic, true, isStochastic, false, p);
                    model.ElementsDictionary[i * faceElements + 78].ElementType = GetSolidElementType(isElastic, true, isStochastic, false, p);
                }
            }
            return elementsPerDirection;
        }

        private static int[] BuildEmbeddedElements(int[] nodesPerDirection)
        {
            int[] elementsPerDirection = new int[] 
                { 
                    (int)((nodesPerDirection[0] - 1) / (nodesPerElementSide[0] - 1)),
                    (int)((nodesPerDirection[1] - 1) / (nodesPerElementSide[1] - 1)),
                    (int)((nodesPerDirection[2] - 1) / (nodesPerElementSide[2] - 1))
                };
            int elements = 0;
            int nX = nodesPerDirection[0];
            int nY = nodesPerDirection[1];
            int nZ = nodesPerDirection[2];

            for (int i = 0; i < elementsPerDirection[2]; i++)
                for (int j = 0; j < elementsPerDirection[1]; j++)
                    for (int k = 0; k < elementsPerDirection[0]; k++)
                    {
                        elements++;
                        Element element = new Element();
                        element.ElementType = GetSolidElementType(false, false, false, true, null);
                        element.ID = elements;
                        element.AddNode(model.NodesDictionary[nX*nY*i + nX*j + k+1]);
                        element.AddNode(model.NodesDictionary[nX*nY*i + nX*j + k+2]);
                        element.AddNode(model.NodesDictionary[nX*nY*i + nX*(j+1) + k+2]);
                        element.AddNode(model.NodesDictionary[nX*nY*i + nX*(j+1) + k+1]);
                        element.AddNode(model.NodesDictionary[nX*nY*(i+1) + nX*j + k+1]);
                        element.AddNode(model.NodesDictionary[nX*nY*(i+1) + nX*j + k+2]);
                        element.AddNode(model.NodesDictionary[nX*nY*(i+1) + nX*(j+1) + k+2]);
                        element.AddNode(model.NodesDictionary[nX*nY*(i+1) + nX*(j+1) + k+1]);
                        model.ElementsDictionary.Add(element.ID, element);
                    }
            return elementsPerDirection;
        }

        private static void BuildSubdomains(int[] elementsPerSubdomainSide, int[] elementsPerDirection, bool hasDampers)
        {
            int[] subdomainsPerDirection = new int[] 
                { 
                    (int)(elementsPerDirection[0] / elementsPerSubdomainSide[0]),
                    (int)(elementsPerDirection[1] / elementsPerSubdomainSide[1]),
                    (int)(elementsPerDirection[2] / elementsPerSubdomainSide[2])
                };
            int subdomains = 0;
            int nX = elementsPerDirection[0];
            int nY = elementsPerDirection[1];
            int nZ = elementsPerDirection[2];

            for (int i = 0; i < subdomainsPerDirection[2]; i++)
                for (int j = 0; j < subdomainsPerDirection[1]; j++)
                    for (int k = 0; k < subdomainsPerDirection[0]; k++)
                    {
                        subdomains++;
                        Subdomain subdomain = new Subdomain() { ID = subdomains };
                        for (int l = 0; l < elementsPerSubdomainSide[2]; l++)
                            for (int m = 0; m < elementsPerSubdomainSide[1]; m++)
                                for (int n = 0; n < elementsPerSubdomainSide[0]; n++)
                                {
                                    int elementID = i*nX*nY*elementsPerSubdomainSide[2] + 
                                        j*nX*elementsPerSubdomainSide[1] +
                                        k*elementsPerSubdomainSide[0] + l*nX*nY + m*nX + n+1;
                                    subdomain.ElementsDictionary.Add(elementID, model.ElementsDictionary[elementID]);
                                }
                        model.SubdomainsDictionary.Add(subdomains, subdomain);
                    }

            if (hasDampers && subdomains > 1) throw new InvalidOperationException("Subdomains and dampers not supported.");
            foreach (var e in model.ElementsDictionary.Where(x => x.Key < 0))
                model.SubdomainsDictionary[1].ElementsDictionary.Add(e.Key, e.Value);
        }

        public static void BuildEmbeddedMesh(int[] nodesPerDirection, double[] distancePerDirection,
            int[] elementsPerSubdomainSide, Surface pinSurface, Surface poreSurface, Surface loadSurface,
            Surface massAccSurface, DOFType loadDirection, double loadAmount, MassAccelerationHistoryLoad massLoad)
        {
            BuildNodes(nodesPerDirection, distancePerDirection, false);
            BuildConstraintsAndLoads(nodesPerDirection, pinSurface, poreSurface, loadSurface,
                loadDirection, loadAmount, false);
            int[] elementsPerDirection = BuildEmbeddedElements(nodesPerDirection);
            BuildElementMassAccelerationLoads(elementsPerDirection, massAccSurface, massLoad);
            BuildSubdomains(elementsPerSubdomainSide, elementsPerDirection, false);
        }

        public static void BuildMesh(int[] nodesPerDirection, double[] distancePerDirection,
            int[] elementsPerSubdomainSide, Surface pinSurface, Surface poreSurface, Surface loadSurface, 
            Surface massAccSurface, DOFType loadDirection, double loadAmount, MassAccelerationHistoryLoad massLoad,
            bool isPorous, bool isElastic, bool isIllConditioned = false, bool isStochastic = false,
            IStochasticMaterialCoefficientsProvider p = null, bool hasDampers = false)
        {
            BuildNodes(nodesPerDirection, distancePerDirection, hasDampers);
            BuildConstraintsAndLoads(nodesPerDirection, pinSurface, poreSurface, loadSurface, loadDirection, loadAmount, hasDampers);
            int[] elementsPerDirection = BuildElements(nodesPerDirection, isPorous, isElastic, isStochastic, isIllConditioned, p, hasDampers);
            BuildElementMassAccelerationLoads(elementsPerDirection, massAccSurface, massLoad); 
            BuildSubdomains(elementsPerSubdomainSide, elementsPerDirection, true);
        }
    }
}
