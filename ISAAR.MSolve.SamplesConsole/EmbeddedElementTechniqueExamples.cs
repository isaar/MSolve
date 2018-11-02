using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ISAAR.MSolve.SamplesConsole
{
    public class EmbeddedElementTechniqueExamples
    {
        public static void SolveEmbeddedElementTechniqueExamples()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Choose model
            EmbeddedExamplesBuilder.ExampleWithEmbedded(model);
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);

            // Skyline Solver
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int totalDOFs = model.TotalDOFs;
            int increments = 100;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[12][DOFType.X],
            model.NodalDOFsDictionary[12][DOFType.Y],
            model.NodalDOFsDictionary[12][DOFType.Z]});

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //Console.WriteLine("checkPoint1 reached");
            Console.WriteLine("Writing results for node 5");
            Console.WriteLine("Dof and Values for Displacement X, Y, Z");
            Console.WriteLine(childAnalyzer.Logs[1][0]);
        }
    }

    public static class EmbeddedExamplesBuilder
    {
        public static void HexaElementsOnly(Model model)
        {
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
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = DOFType.X });
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = DOFType.Y });
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = DOFType.Z });
                    nodeID++;
                }
            }

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };

            // first element definition
            Element e1;
            int ID2 = 1;
            e1 = new Element()
            {
                ID = 1,
                ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
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

            // second element definition
            Element e2 = new Element()
            {
                ID = 2,
                ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
            };
            ID2 = 1;
            for (int j = 0; j < 2; j++)
            {
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 5, model.NodesDictionary[4 * (ID2 - 1) + 5]);
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 6, model.NodesDictionary[4 * (ID2 - 1) + 6]);
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 8, model.NodesDictionary[4 * (ID2 - 1) + 8]);
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 7, model.NodesDictionary[4 * (ID2 - 1) + 7]);
                ID2++;
            }
            model.ElementsDictionary.Add(e2.ID, e2);

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

        public static void BeamElementOnly(Model model)
        {
            // define mechanical properties properties
            double youngModulus = 16710.0;
            double poissonRatio = 0.034;
            double nodalLoad = 50.0;
            double area = 5.594673861218848e-003;
            double inertiaY = 2.490804749753243e-006;
            double inertiaZ = 2.490804749753243e-006;
            double torsionalInertia = inertiaY / 2.0;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;

            // geometry
            model.NodesDictionary.Add(13, new Node() { ID = 13, X = 0.5, Y = 0.5, Z = 0.5 });
            model.NodesDictionary.Add(14, new Node() { ID = 14, X = 0.5, Y = 0.5, Z = 1.5 });

            // constraints
            model.NodesDictionary[13].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[13].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[13].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[13].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[13].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[13].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            model.NodesDictionary[14].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[14].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[14].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            // Create new 3D material
            ElasticMaterial3D material_1 = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Create new Beam3D section and element
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
            // element nodes
            IList<Node> elementNodes = new List<Node>();
            elementNodes.Add(model.NodesDictionary[13]);
            elementNodes.Add(model.NodesDictionary[14]);
            var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, material_1, 7.85, beamSection);
            var element = new Element { ID = 3, ElementType = beam_1 };
            int subdomainID = 1;
            
            element.NodesDictionary.Add(13, model.NodesDictionary[13]);
            element.NodesDictionary.Add(14, model.NodesDictionary[14]);

            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(element.ID, element);

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[14], DOF = DOFType.Z });
        }

        public static void ExampleWithEmbedded(Model model)
        {
            HexaElementsOnly(model);
            BeamElementOnly(model);
            model.Loads.RemoveAt(4);
            model.NodesDictionary[13].Constraints.Remove(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[13].Constraints.Remove(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[13].Constraints.Remove(new Constraint { DOF = DOFType.Z });
            var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key < 3).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key >= 3).Select(kv => kv.Value));
        }
    }
}
