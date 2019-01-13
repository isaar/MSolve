using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Mesh.SuiteSparse;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests.Visualization
{
    class VTKWriterTests
    {

        public static void Main()
        {
            Model2D model = CreateModel();
            var writer = new VTKWriter(model);
            writer.InitializeFile("test1", false);
            writer.WriteScalarField("Distance_from_orgin", GenerateScalarField(model));
            writer.CloseCurrentFile();
        }

        private static Model2D CreateModel()
        {
            // Mesh
            Model2D model = new Model2D();
            var meshGenerator = new UniformRectilinearMeshGenerator(100, 50, 40, 25);
            (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) = meshGenerator.CreateMesh();

            // Nodes
            foreach (XNode2D node in nodes) model.AddNode(node);

            // Elements
            var integration = new XSimpleIntegration2D();
            //var integration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //        new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            //var jIntegration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order4x4,
            //    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4));

            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(2e6, 0.3);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration));
            }

            return model;
        }

        private static double[] GenerateScalarField(Model2D model)
        {
            double[] nodalValues = new double[model.Nodes.Count];
            for (int n = 0; n < model.Nodes.Count; ++ n)
            {
                XNode2D node = model.Nodes[n];
                nodalValues[n] = Math.Sqrt(node.X * node.X + node.Y * node.Y);
            }
            return nodalValues;
        }
    }
}
