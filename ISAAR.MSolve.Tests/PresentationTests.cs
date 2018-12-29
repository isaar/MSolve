using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.GMSH;

namespace ISAAR.MSolve.Tests
{
    public class PresentationTests
    {
        //[Fact]
        private static void SolveStaticQuadRetainingWall()
        {
            #region Model
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 2.1e09;
            double poissonRatio = 0.3;

            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            Model model = new Model();

            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            #region Nodes
            model.NodesDictionary.Add(0, new Node { ID = 0, X = 0.0, Y = 0.0, Z = 0.0 });
            model.NodesDictionary.Add(1, new Node { ID = 1, X = 0.6, Y = 0.0, Z = 0.0 });

            model.NodesDictionary.Add(2, new Node { ID = 2, X = 0.0, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(3, new Node { ID = 3, X = 0.6, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(4, new Node { ID = 4, X = 1.2, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(5, new Node { ID = 5, X = 1.8, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(6, new Node { ID = 6, X = 2.4, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(7, new Node { ID = 7, X = 3.0, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(8, new Node { ID = 8, X = 3.6, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(9, new Node { ID = 9, X = 4.3, Y = 0.5, Z = 0.0 });

            model.NodesDictionary.Add(10, new Node { ID = 10, X = 0.0, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(11, new Node { ID = 11, X = 0.6, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(12, new Node { ID = 12, X = 1.2, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(13, new Node { ID = 13, X = 1.8, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(14, new Node { ID = 14, X = 2.4, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(15, new Node { ID = 15, X = 3.0, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(16, new Node { ID = 16, X = 3.6, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(17, new Node { ID = 17, X = 4.3, Y = 1.2, Z = 0.0 });

            model.NodesDictionary.Add(18, new Node { ID = 18, X = 3, Y = 1.9, Z = 0.0 });
            model.NodesDictionary.Add(19, new Node { ID = 19, X = 3.5756, Y = 1.9, Z = 0.0 });

            model.NodesDictionary.Add(20, new Node { ID = 20, X = 3, Y = 2.6, Z = 0.0 });
            model.NodesDictionary.Add(21, new Node { ID = 21, X = 3.5512, Y = 2.6, Z = 0.0 });

            model.NodesDictionary.Add(22, new Node { ID = 22, X = 3, Y = 3.3, Z = 0.0 });
            model.NodesDictionary.Add(23, new Node { ID = 23, X = 3.5267, Y = 3.3, Z = 0.0 });

            model.NodesDictionary.Add(24, new Node { ID = 24, X = 3, Y = 4.0, Z = 0.0 });
            model.NodesDictionary.Add(25, new Node { ID = 25, X = 3.5023, Y = 4.0, Z = 0.0 });

            model.NodesDictionary.Add(26, new Node { ID = 26, X = 3, Y = 4.7, Z = 0.0 });
            model.NodesDictionary.Add(27, new Node { ID = 27, X = 3.4779, Y = 4.7, Z = 0.0 });

            model.NodesDictionary.Add(28, new Node { ID = 28, X = 3, Y = 5.5, Z = 0.0 });
            model.NodesDictionary.Add(29, new Node { ID = 29, X = 3.45, Y = 5.5, Z = 0.0 });
            #endregion

            #region Elements

            #region element0
            var element0 = new Element()
            {
                ID = 0,
                ElementType = new Quad4(material)
            };
            element0.AddNode(model.NodesDictionary[0]);
            element0.AddNode(model.NodesDictionary[1]);
            element0.AddNode(model.NodesDictionary[3]);
            element0.AddNode(model.NodesDictionary[2]);
            model.ElementsDictionary.Add(element0.ID, element0);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element0.ID, element0);
            #endregion

            #region element1
            var element1 = new Element()
            {
                ID = 1,
                ElementType = new Quad4(material)
            };
            element1.AddNode(model.NodesDictionary[2]);
            element1.AddNode(model.NodesDictionary[3]);
            element1.AddNode(model.NodesDictionary[11]);
            element1.AddNode(model.NodesDictionary[10]);
            model.ElementsDictionary.Add(element1.ID, element1);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element1.ID, element1);
            #endregion

            #region element2
            var element2 = new Element()
            {
                ID = 2,
                ElementType = new Quad4(material)
            };
            element2.AddNode(model.NodesDictionary[3]);
            element2.AddNode(model.NodesDictionary[4]);
            element2.AddNode(model.NodesDictionary[12]);
            element2.AddNode(model.NodesDictionary[11]);
            model.ElementsDictionary.Add(element2.ID, element2);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element2.ID, element2);
            #endregion

            #region element3
            var element3 = new Element()
            {
                ID = 3,
                ElementType = new Quad4(material)
            };
            element3.AddNode(model.NodesDictionary[4]);
            element3.AddNode(model.NodesDictionary[5]);
            element3.AddNode(model.NodesDictionary[13]);
            element3.AddNode(model.NodesDictionary[12]);
            model.ElementsDictionary.Add(element3.ID, element3);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element3.ID, element3);
            #endregion

            #region element4
            var element4 = new Element()
            {
                ID = 4,
                ElementType = new Quad4(material)
            };
            element4.AddNode(model.NodesDictionary[5]);
            element4.AddNode(model.NodesDictionary[6]);
            element4.AddNode(model.NodesDictionary[14]);
            element4.AddNode(model.NodesDictionary[13]);
            model.ElementsDictionary.Add(element4.ID, element4);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element4.ID, element4);
            #endregion

            #region element5
            var element5 = new Element()
            {
                ID = 5,
                ElementType = new Quad4(material)
            };
            element5.AddNode(model.NodesDictionary[6]);
            element5.AddNode(model.NodesDictionary[7]);
            element5.AddNode(model.NodesDictionary[15]);
            element5.AddNode(model.NodesDictionary[14]);
            model.ElementsDictionary.Add(element5.ID, element5);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element5.ID, element5);
            #endregion

            #region element6
            var element6 = new Element()
            {
                ID = 6,
                ElementType = new Quad4(material)
            };
            element6.AddNode(model.NodesDictionary[7]);
            element6.AddNode(model.NodesDictionary[8]);
            element6.AddNode(model.NodesDictionary[16]);
            element6.AddNode(model.NodesDictionary[15]);
            model.ElementsDictionary.Add(element6.ID, element6);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element6.ID, element6);
            #endregion

            #region element7
            var element7 = new Element()
            {
                ID = 7,
                ElementType = new Quad4(material)
            };
            element7.AddNode(model.NodesDictionary[8]);
            element7.AddNode(model.NodesDictionary[9]);
            element7.AddNode(model.NodesDictionary[17]);
            element7.AddNode(model.NodesDictionary[16]);
            model.ElementsDictionary.Add(element7.ID, element7);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element7.ID, element7);
            #endregion

            #region element8
            var element8 = new Element()
            {
                ID = 8,
                ElementType = new Quad4(material)
            };
            element8.AddNode(model.NodesDictionary[15]);
            element8.AddNode(model.NodesDictionary[16]);
            element8.AddNode(model.NodesDictionary[19]);
            element8.AddNode(model.NodesDictionary[18]);
            model.ElementsDictionary.Add(element8.ID, element8);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element8.ID, element8);
            #endregion

            #region element9
            var element9 = new Element()
            {
                ID = 9,
                ElementType = new Quad4(material)
            };
            element9.AddNode(model.NodesDictionary[18]);
            element9.AddNode(model.NodesDictionary[19]);
            element9.AddNode(model.NodesDictionary[21]);
            element9.AddNode(model.NodesDictionary[20]);
            model.ElementsDictionary.Add(element9.ID, element9);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element9.ID, element9);
            #endregion

            #region element10
            var element10 = new Element()
            {
                ID = 10,
                ElementType = new Quad4(material)
            };
            element10.AddNode(model.NodesDictionary[20]);
            element10.AddNode(model.NodesDictionary[21]);
            element10.AddNode(model.NodesDictionary[23]);
            element10.AddNode(model.NodesDictionary[22]);
            model.ElementsDictionary.Add(element10.ID, element10);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element10.ID, element10);
            #endregion

            #region element11
            var element11 = new Element()
            {
                ID = 11,
                ElementType = new Quad4(material)
            };
            element11.AddNode(model.NodesDictionary[22]);
            element11.AddNode(model.NodesDictionary[23]);
            element11.AddNode(model.NodesDictionary[25]);
            element11.AddNode(model.NodesDictionary[24]);
            model.ElementsDictionary.Add(element11.ID, element11);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element11.ID, element11);
            #endregion

            #region element12
            var element12 = new Element()
            {
                ID = 12,
                ElementType = new Quad4(material)
            };
            element12.AddNode(model.NodesDictionary[24]);
            element12.AddNode(model.NodesDictionary[25]);
            element12.AddNode(model.NodesDictionary[27]);
            element12.AddNode(model.NodesDictionary[26]);
            model.ElementsDictionary.Add(element12.ID, element12);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element12.ID, element12);
            #endregion

            #region element13
            var element13 = new Element()
            {
                ID = 13,
                ElementType = new Quad4(material)
            };
            element13.AddNode(model.NodesDictionary[26]);
            element13.AddNode(model.NodesDictionary[27]);
            element13.AddNode(model.NodesDictionary[29]);
            element13.AddNode(model.NodesDictionary[28]);
            model.ElementsDictionary.Add(element13.ID, element13);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element13.ID, element13);
            #endregion
            #endregion

            #region Constrains
            var constrainedNodes = new int[] { 0, 1, 3, 4, 5, 6, 7, 8, 9 };
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = DOFType.Y });
            }

            #endregion

            #region Loads
            //ground vertical loads
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[10], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[11], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[12], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[13], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[14], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[15], DOF = DOFType.Y });

            //ground horizontal loads
            model.Loads.Add(new Load() { Amount = -2130, Node = model.NodesDictionary[28], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -11490, Node = model.NodesDictionary[26], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -20990, Node = model.NodesDictionary[24], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -30790, Node = model.NodesDictionary[22], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -40600, Node = model.NodesDictionary[22], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[20], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -50390, Node = model.NodesDictionary[18], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -28460, Node = model.NodesDictionary[15], DOF = DOFType.X });
            #endregion

            model.ConnectDataStructures();

            #endregion

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

        }

        //[Fact]
        private static void SolveDynamicQuadRetainingWall()
        {
            #region Model
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 2.1e09;
            double poissonRatio = 0.3;
            double density = 20;

            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            Model model = new Model();

            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            #region Nodes
            model.NodesDictionary.Add(0, new Node { ID = 0, X = 0.0, Y = 0.0, Z = 0.0 });
            model.NodesDictionary.Add(1, new Node { ID = 1, X = 0.6, Y = 0.0, Z = 0.0 });

            model.NodesDictionary.Add(2, new Node { ID = 2, X = 0.0, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(3, new Node { ID = 3, X = 0.6, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(4, new Node { ID = 4, X = 1.2, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(5, new Node { ID = 5, X = 1.8, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(6, new Node { ID = 6, X = 2.4, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(7, new Node { ID = 7, X = 3.0, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(8, new Node { ID = 8, X = 3.6, Y = 0.5, Z = 0.0 });
            model.NodesDictionary.Add(9, new Node { ID = 9, X = 4.3, Y = 0.5, Z = 0.0 });

            model.NodesDictionary.Add(10, new Node { ID = 10, X = 0.0, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(11, new Node { ID = 11, X = 0.6, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(12, new Node { ID = 12, X = 1.2, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(13, new Node { ID = 13, X = 1.8, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(14, new Node { ID = 14, X = 2.4, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(15, new Node { ID = 15, X = 3.0, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(16, new Node { ID = 16, X = 3.6, Y = 1.2, Z = 0.0 });
            model.NodesDictionary.Add(17, new Node { ID = 17, X = 4.3, Y = 1.2, Z = 0.0 });

            model.NodesDictionary.Add(18, new Node { ID = 18, X = 3, Y = 1.9, Z = 0.0 });
            model.NodesDictionary.Add(19, new Node { ID = 19, X = 3.5756, Y = 1.9, Z = 0.0 });

            model.NodesDictionary.Add(20, new Node { ID = 20, X = 3, Y = 2.6, Z = 0.0 });
            model.NodesDictionary.Add(21, new Node { ID = 21, X = 3.5512, Y = 2.6, Z = 0.0 });

            model.NodesDictionary.Add(22, new Node { ID = 22, X = 3, Y = 3.3, Z = 0.0 });
            model.NodesDictionary.Add(23, new Node { ID = 23, X = 3.5267, Y = 3.3, Z = 0.0 });

            model.NodesDictionary.Add(24, new Node { ID = 24, X = 3, Y = 4.0, Z = 0.0 });
            model.NodesDictionary.Add(25, new Node { ID = 25, X = 3.5023, Y = 4.0, Z = 0.0 });

            model.NodesDictionary.Add(26, new Node { ID = 26, X = 3, Y = 4.7, Z = 0.0 });
            model.NodesDictionary.Add(27, new Node { ID = 27, X = 3.4779, Y = 4.7, Z = 0.0 });

            model.NodesDictionary.Add(28, new Node { ID = 28, X = 3, Y = 5.5, Z = 0.0 });
            model.NodesDictionary.Add(29, new Node { ID = 29, X = 3.45, Y = 5.5, Z = 0.0 });
            #endregion

            #region Elements

            #region element0
            var element0 = new Element()
            {
                ID = 0,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element0.AddNode(model.NodesDictionary[0]);
            element0.AddNode(model.NodesDictionary[1]);
            element0.AddNode(model.NodesDictionary[3]);
            element0.AddNode(model.NodesDictionary[2]);
            model.ElementsDictionary.Add(element0.ID, element0);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element0.ID, element0);
            #endregion

            #region element1
            var element1 = new Element()
            {
                ID = 1,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element1.AddNode(model.NodesDictionary[2]);
            element1.AddNode(model.NodesDictionary[3]);
            element1.AddNode(model.NodesDictionary[11]);
            element1.AddNode(model.NodesDictionary[10]);
            model.ElementsDictionary.Add(element1.ID, element1);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element1.ID, element1);
            #endregion

            #region element2
            var element2 = new Element()
            {
                ID = 2,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element2.AddNode(model.NodesDictionary[3]);
            element2.AddNode(model.NodesDictionary[4]);
            element2.AddNode(model.NodesDictionary[12]);
            element2.AddNode(model.NodesDictionary[11]);
            model.ElementsDictionary.Add(element2.ID, element2);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element2.ID, element2);
            #endregion

            #region element3
            var element3 = new Element()
            {
                ID = 3,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element3.AddNode(model.NodesDictionary[4]);
            element3.AddNode(model.NodesDictionary[5]);
            element3.AddNode(model.NodesDictionary[13]);
            element3.AddNode(model.NodesDictionary[12]);
            model.ElementsDictionary.Add(element3.ID, element3);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element3.ID, element3);
            #endregion

            #region element4
            var element4 = new Element()
            {
                ID = 4,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element4.AddNode(model.NodesDictionary[5]);
            element4.AddNode(model.NodesDictionary[6]);
            element4.AddNode(model.NodesDictionary[14]);
            element4.AddNode(model.NodesDictionary[13]);
            model.ElementsDictionary.Add(element4.ID, element4);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element4.ID, element4);
            #endregion

            #region element5
            var element5 = new Element()
            {
                ID = 5,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element5.AddNode(model.NodesDictionary[6]);
            element5.AddNode(model.NodesDictionary[7]);
            element5.AddNode(model.NodesDictionary[15]);
            element5.AddNode(model.NodesDictionary[14]);
            model.ElementsDictionary.Add(element5.ID, element5);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element5.ID, element5);
            #endregion

            #region element6
            var element6 = new Element()
            {
                ID = 6,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element6.AddNode(model.NodesDictionary[7]);
            element6.AddNode(model.NodesDictionary[8]);
            element6.AddNode(model.NodesDictionary[16]);
            element6.AddNode(model.NodesDictionary[15]);
            model.ElementsDictionary.Add(element6.ID, element6);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element6.ID, element6);
            #endregion

            #region element7
            var element7 = new Element()
            {
                ID = 7,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element7.AddNode(model.NodesDictionary[8]);
            element7.AddNode(model.NodesDictionary[9]);
            element7.AddNode(model.NodesDictionary[17]);
            element7.AddNode(model.NodesDictionary[16]);
            model.ElementsDictionary.Add(element7.ID, element7);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element7.ID, element7);
            #endregion

            #region element8
            var element8 = new Element()
            {
                ID = 8,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element8.AddNode(model.NodesDictionary[15]);
            element8.AddNode(model.NodesDictionary[16]);
            element8.AddNode(model.NodesDictionary[19]);
            element8.AddNode(model.NodesDictionary[18]);
            model.ElementsDictionary.Add(element8.ID, element8);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element8.ID, element8);
            #endregion

            #region element9
            var element9 = new Element()
            {
                ID = 9,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element9.AddNode(model.NodesDictionary[18]);
            element9.AddNode(model.NodesDictionary[19]);
            element9.AddNode(model.NodesDictionary[21]);
            element9.AddNode(model.NodesDictionary[20]);
            model.ElementsDictionary.Add(element9.ID, element9);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element9.ID, element9);
            #endregion

            #region element10
            var element10 = new Element()
            {
                ID = 10,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element10.AddNode(model.NodesDictionary[20]);
            element10.AddNode(model.NodesDictionary[21]);
            element10.AddNode(model.NodesDictionary[23]);
            element10.AddNode(model.NodesDictionary[22]);
            model.ElementsDictionary.Add(element10.ID, element10);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element10.ID, element10);
            #endregion

            #region element11
            var element11 = new Element()
            {
                ID = 11,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element11.AddNode(model.NodesDictionary[22]);
            element11.AddNode(model.NodesDictionary[23]);
            element11.AddNode(model.NodesDictionary[25]);
            element11.AddNode(model.NodesDictionary[24]);
            model.ElementsDictionary.Add(element11.ID, element11);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element11.ID, element11);
            #endregion

            #region element12
            var element12 = new Element()
            {
                ID = 12,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element12.AddNode(model.NodesDictionary[24]);
            element12.AddNode(model.NodesDictionary[25]);
            element12.AddNode(model.NodesDictionary[27]);
            element12.AddNode(model.NodesDictionary[26]);
            model.ElementsDictionary.Add(element12.ID, element12);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element12.ID, element12);
            #endregion

            #region element13
            var element13 = new Element()
            {
                ID = 13,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element13.AddNode(model.NodesDictionary[26]);
            element13.AddNode(model.NodesDictionary[27]);
            element13.AddNode(model.NodesDictionary[29]);
            element13.AddNode(model.NodesDictionary[28]);
            model.ElementsDictionary.Add(element13.ID, element13);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element13.ID, element13);
            #endregion
            #endregion

            #region Constrains
            var constrainedNodes = new int[] { 0, 1, 3, 4, 5, 6, 7, 8, 9 };
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = DOFType.Y });
            }

            #endregion

            #region Loads
            //ground vertical loads
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[10], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[11], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[12], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[13], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[14], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[15], DOF = DOFType.Y });

            //ground horizontal loads
            model.Loads.Add(new Load() { Amount = -2130, Node = model.NodesDictionary[28], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -11490, Node = model.NodesDictionary[26], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -20990, Node = model.NodesDictionary[24], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -30790, Node = model.NodesDictionary[22], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -40600, Node = model.NodesDictionary[22], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[20], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -50390, Node = model.NodesDictionary[18], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = -28460, Node = model.NodesDictionary[15], DOF = DOFType.X });

            model.MassAccelerationHistoryLoads.Add(new MassAccelerationHistoryLoad("..\\..\\..\\elcentro_NS.dat", 1) { DOF = DOFType.X });
            #endregion

            model.ConnectDataStructures();

            #endregion

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.6, 1, 0.02, 53.74);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

        }

        //[Fact]
        private static void SolveStaticLinearWall()
        {
            #region Read Data
            string workingDirectory = @"C:\Users\Dimitris\Desktop\Presentation";
            string meshFileName = "wall.msh";

            (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) = GenerateMeshFromGmsh(workingDirectory, meshFileName);
            #endregion

            #region CreateModel

            const double height = 3.5;
            const double thickness = 0.1;
            const double youngModulus = 2E6;
            const double poissonRatio = 0.3;
            const double maxLoad = 1000.0;
            // Initialize
            int numberOfNodes = nodes.Count;
            int numberOfElements = elements.Count;
            VectorExtensions.AssignTotalAffinityCount();

            // Materials
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };

            // Subdomains
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Nodes
            for (int i = 0; i < numberOfNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            for (int i = 0; i < numberOfElements; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, elementWrapper);
            }

            // Constraints
            double tol = 1E-10;
            Node2D[] constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol).ToArray();
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = DOFType.X });
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = DOFType.Y });
            }

            // Loads
            Node2D[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.Loads.Add(new Load() { Amount = maxLoad, Node = loadedNodes[0], DOF = DOFType.X });

            // Finalize
            model.ConnectDataStructures();


            #endregion

            #region Analysis
            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            // Logging displacement, strain, and stress fields.
            string outputDirectory = workingDirectory + "\\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion
        }

        //[Fact]
        private static void SolveStaticNonLinearWall()
        {
            #region Read Data
            string workingDirectory = @"C:\Users\Dimitris\Desktop\Presentation";
            string meshFileName = "wall.msh";

            (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) = GenerateMeshFromGmsh(workingDirectory, meshFileName);
            #endregion

            #region CreateModel

            const double height = 3.5;
            const double thickness = 0.1;
            const double youngModulus = 2E6;
            const double poissonRatio = 0.3;
            const double maxLoad = 1000.0;
            // Initialize
            int numberOfNodes = nodes.Count;
            int numberOfElements = elements.Count;
            VectorExtensions.AssignTotalAffinityCount();

            // Materials
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };



            // Subdomains
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Nodes
            for (int i = 0; i < numberOfNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            for (int i = 0; i < numberOfElements; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, elementWrapper);
            }

            // Constraints
            double tol = 1E-10;
            Node2D[] constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol).ToArray();
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = DOFType.X });
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = DOFType.Y });
            }

            // Loads
            Node2D[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.Loads.Add(new Load() { Amount = maxLoad, Node = loadedNodes[0], DOF = DOFType.X });

            // Finalize
            model.ConnectDataStructures();


            #endregion

            #region Analysis
            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            // Logging displacement, strain, and stress fields.
            string outputDirectory = workingDirectory + "\\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion
        }

        //[Fact]
        private static void SolveDynamicWall()
        {
            #region Read Data
            string workingDirectory = @"C:\Users\Dimitris\Desktop\Presentation";
            string meshFileName = "wall.msh";

            (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) = GenerateMeshFromGmsh(workingDirectory, meshFileName);
            #endregion

            #region CreateModel

            const double height = 3.5;
            const double thickness = 0.1;
            const double youngModulus = 2E6;
            const double poissonRatio = 0.3;
            const double maxLoad = 1000.0;
            // Initialize
            int numberOfNodes = nodes.Count;
            int numberOfElements = elements.Count;
            VectorExtensions.AssignTotalAffinityCount();

            // Materials
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };

            DynamicMaterial dynamicMaterial = new DynamicMaterial(25, 0.05, 0.05);

            // Subdomains
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Nodes
            for (int i = 0; i < numberOfNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, dynamicMaterial);
            for (int i = 0; i < numberOfElements; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, elementWrapper);
            }

            // Constraints
            double tol = 1E-10;
            Node2D[] constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol).ToArray();
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = DOFType.X });
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = DOFType.Y });
            }

            // Loads
            Node2D[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.Loads.Add(new Load() { Amount = maxLoad, Node = loadedNodes[0], DOF = DOFType.X });

            model.MassAccelerationHistoryLoads.Add(new MassAccelerationHistoryLoad("..\\..\\..\\elcentro_NS.dat", 1) { DOF = DOFType.X });


            // Finalize
            model.ConnectDataStructures();


            #endregion

            #region Analysis
            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.6, 1, 0.02, 53.74);


            // Logging displacement, strain, and stress fields.
            string outputDirectory = workingDirectory + "\\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion
        }

        private static (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) GenerateMeshFromGmsh(string workingDirectory, string meshFileName)
        {
            using (var reader = new GmshReader2D(workingDirectory + "\\" + meshFileName))
            {
                return reader.CreateMesh();
            }
        }


        //[Fact]
        public void Solve4QuadProblem()
        {
            VectorExtensions.AssignTotalAffinityCount();

            #region Nodes
            var node0 = new Node { ID = 0, X = 0.0, Y = 0.0 };
            var node1 = new Node { ID = 1, X = 2.0, Y = 0.0 };
            var node2 = new Node { ID = 2, X = 4.0, Y = 0.0 };
            var node3 = new Node { ID = 3, X = 0.0, Y = 2.0 };
            var node4 = new Node { ID = 4, X = 2.0, Y = 2.0 };
            var node5 = new Node { ID = 5, X = 4.0, Y = 2.0 };
            var node6 = new Node { ID = 6, X = 0.0, Y = 4.0 };
            var node7 = new Node { ID = 7, X = 2.0, Y = 4.0 };
            var node8 = new Node { ID = 8, X = 4.0, Y = 4.0 };
            #endregion

            #region material
            var material = new ElasticMaterial2D(StressState2D.PlaneStress);
            material.YoungModulus = 3e7;
            material.PoissonRatio = 0.2;

            double thickness = 0.3;
            #endregion

            #region elements

            var element0 = new Element { ID = 0, ElementType = new Quad4(material) { Thickness = thickness } };
            element0.AddNodes(new[] { node0, node1, node4, node3 });

            var element1 = new Element { ID = 1, ElementType = new Quad4(material) { Thickness = thickness } };
            element1.AddNodes(new[] { node1, node2, node5, node4 });

            var element2 = new Element { ID = 2, ElementType = new Quad4(material) { Thickness = thickness } };
            element2.AddNodes(new[] { node3, node4, node7, node6 });

            var element3 = new Element { ID = 3, ElementType = new Quad4(material) { Thickness = thickness } };
            element3.AddNodes(new[] { node4, node5, node8, node7 });
            #endregion

            #region loads
            var load0 = new Load { Amount = 50, DOF = DOFType.X, Node = node8 };
            #endregion

            #region constraints
            node0.Constraints.Add(new Constraint { DOF = DOFType.X });
            node0.Constraints.Add(new Constraint { DOF = DOFType.Y });

            node2.Constraints.Add(new Constraint { DOF = DOFType.X });
            node2.Constraints.Add(new Constraint { DOF = DOFType.Y });
            #endregion

            #region subdomains
            var subdomain0 = new Subdomain { ID = 0 };
            subdomain0.ElementsDictionary.Add(0, element0);
            subdomain0.ElementsDictionary.Add(1, element1);
            subdomain0.ElementsDictionary.Add(2, element2);
            subdomain0.ElementsDictionary.Add(3, element3);
            #endregion

            #region model
            var model = new Model();

            model.NodesDictionary.Add(0, node0);
            model.NodesDictionary.Add(1, node1);
            model.NodesDictionary.Add(2, node2);
            model.NodesDictionary.Add(3, node3);
            model.NodesDictionary.Add(4, node4);
            model.NodesDictionary.Add(5, node5);
            model.NodesDictionary.Add(6, node6);
            model.NodesDictionary.Add(7, node7);
            model.NodesDictionary.Add(8, node8);

            model.ElementsDictionary.Add(0, element0);
            model.ElementsDictionary.Add(1, element1);
            model.ElementsDictionary.Add(2, element2);
            model.ElementsDictionary.Add(3, element3);

            model.Loads.Add(load0);

            model.SubdomainsDictionary.Add(0, subdomain0);

            model.ConnectDataStructures();
            #endregion

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            var solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            var childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            var parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
