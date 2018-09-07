using System;
using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Gmsh
{
    class GmshReaderCantilever
    {
        private static readonly string directory =
           Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Resources\\";
        private static readonly string filename = "gmsh_cantilevel_test.msh";
        private static readonly double E = 210e6; // kN/m^2
        private static readonly double v = 0.3;
        private static readonly double length = 5.0; //m
        private static readonly double depth = 0.5; //m
        private static readonly double thickness = 0.5; //m
        private static readonly double I = thickness * Math.Pow(depth, 3.0) / 12.0;
        private static readonly double load = -10; //kN
        private static readonly double analyticDisplacement = load * Math.Pow(length, 3.0) / (3.0 * E * I);

        public static void Main()
        {
            var test = new GmshReaderCantilever();
            Model2D model = test.CreateModel();
            Tuple<XNode2D, double> maxDisplacement = test.FindMaxDisplacement(model);
            Console.WriteLine("Expected max displacement: " + analyticDisplacement);
            Console.WriteLine("Computed max displacement: " + maxDisplacement.Item2 + " at node " + maxDisplacement.Item1);
        }

        private Model2D CreateModel()
        {
            var reader = new GmshReader(filename);
            Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> meshEntites = reader.ReadMesh();
            Model2D model = new Model2D();

            foreach (var node in meshEntites.Item1)
            {
                model.AddNode(node);
            }
            var integration = new XSimpleIntegration2D();
            foreach (var el in meshEntites.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, thickness);
                var element = new XContinuumElement2D(el.ElementType, el.Nodes, materialField, integration);
                model.AddElement(element);
            }

            ApplyBoundaryConditions(model);
            return model;
        }

        private void ApplyBoundaryConditions(Model2D model)
        {
            var finder = new EntityFinder(model, 1e-6);
            IReadOnlyList<XNode2D> leftNodes = finder.FindNodesWithX(0.0);
            XNode2D topRightNode = finder.FindNodeWith(length, depth);

            // Constrain the nodes of the left edge
            foreach (var node in leftNodes)
            {
                model.AddConstraint(node, DisplacementDof.X, 0.0);
                model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }

            // Apply the load on the top right node
            model.AddNodalLoad(topRightNode, DisplacementDof.Y, load);
        }

        public Tuple<XNode2D, double> FindMaxDisplacement(Model2D model)
        {
            var solver = new SkylineSolver(model);
            solver.Initialize();
            solver.Solve();

            double min = double.MaxValue;
            XNode2D mostDisplaced = null;
            for (int dof = 0; dof < solver.Solution.Length; ++dof)
            {
                if (solver.Solution[dof] < min)
                {
                    min = solver.Solution[dof];
                    int nodeID = dof / 2; // node i: dofX=2i, dofY=2i+1
                    mostDisplaced = model.Nodes[nodeID];
                }
            }
            return new Tuple<XNode2D, double>(mostDisplaced, min);
        }

    }
}
