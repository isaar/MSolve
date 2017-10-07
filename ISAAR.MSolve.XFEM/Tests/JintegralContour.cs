using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.XFEM.Analysis;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Output;
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests
{
    class JintegralContour
    {
        public static void Main()
        {
            int elementsPerAxis = 15;
            double width = elementsPerAxis * 2;
            double radius = width / 5;
            double center = width / 2;
            string outputFile = "weights";

            // Create model
            var model = new Model2D();
            var meshGenerator = new UniformRectilinearMeshGenerator(width, width, elementsPerAxis, elementsPerAxis);
            Tuple<XNode2D[], List<XNode2D[]>> meshEntities =  meshGenerator.CreateMesh();
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);
            var integration = new XSimpleIntegration2D();
            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStress(2.1e6, 0.3, 1.0);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration));
            }

            // Assign weights
            var circle = new Circle2D(new CartesianPoint2D(center, center), radius);
            double[] weights = new double[model.Nodes.Count];
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                CirclePointPosition pos = circle.FindRelativePositionOfPoint(model.Nodes[i]);
                if (pos == CirclePointPosition.Inside) weights[i] = 1.0;
            }

            // Print results
            var writer = new VTKWriter(model);
            writer.InitializeFile(outputFile);
            writer.WriteScalarField("weights", weights);
            writer.CloseCurrentFile();
        }
    }
}
