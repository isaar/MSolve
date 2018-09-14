using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class Fillet
    {
        private static readonly double bottomWidth = 375.0, topWidth = 70.0, radius = 20.0;// mm
        private static readonly double flangeHeight = 75.0, totalHeight = 150.0; // mm
        private static readonly double crackHeight = flangeHeight + radius, crackLength = 5; //mm
        private static readonly double E = 2.1e12; // kN/mm^2
        private static readonly double v = 0.3, t = 1.0; // mm
        private static readonly double load = 20.0; // kN
        private static readonly bool fixBottom = false;
        private static readonly string triMeshFile = "fillet_tri.msh";
        private static readonly string quadMeshFile = "fillet_quad.msh";
        private static readonly string outputMesh = "fillet_mesh";
        private static readonly string outputCrack = "fillet_crack.vtk";
        private static readonly bool triangularMesh = false;
        private static readonly bool integrationWithTriangles = false;
        
        public static void Main()
        {
            // Parameters
            double jIntegralRadiusOverElementSize = 1.25;
            double fractureToughness = double.MaxValue;
            int maxIterations = 12;
            double growthLength = 5; // mm

            var benchmark = new Fillet(jIntegralRadiusOverElementSize, fractureToughness,
                maxIterations, growthLength);
            IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze();
            Console.WriteLine("Crack path:");
            foreach (var point in crackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }

            Plot(benchmark.model, crackPath);
        }

        private static void Plot(Model2D model, IReadOnlyList<ICartesianPoint2D> crackPath)
        {
            var writer = new VTKWriter(model);
            writer.InitializeFile(outputMesh, false);
        }

        public Model2D model;
        private TrackingExteriorCrackLSM crack;
        private IIntegrationStrategy2D<XContinuumElement2D> integration;
        private IIntegrationStrategy2D<XContinuumElement2D> jIntegration;
        private readonly double jIntegralRadiusOverElementSize;
        private Propagator propagator;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private readonly double fractureToughness;
        private readonly int maxIterations;
        private readonly double growthLength;

        public Fillet(double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations, double growthLength)
        {
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
            this.growthLength = growthLength;

            CreateModel();
        }

        private void CreateModel()
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, t);
            model = new Model2D();
            HandleIntegrations();
            CreateMesh();
            ApplyBoundaryConditions();
            HandleCrack();
        }

        private void HandleIntegrations()
        {
            if (integrationWithTriangles)
            {
                ITriangulator2D triangulator = new IncrementalTriangulator();
                integration = new IntegrationForCrackPropagation2D(
                    new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order2Points3, crack, triangulator),
                    new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order2Points3, crack, triangulator));
                jIntegration = new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order3Points4, crack,
                    triangulator);
            }
            else
            {
                integration = new IntegrationForCrackPropagation2D(
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
                jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);
                //var jIntegration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
                //    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4));
            }
        }

        private void CreateMesh()
        {
            GmshReader meshReader;
            if (triangularMesh) meshReader = new GmshReader(triMeshFile);
            else meshReader = new GmshReader(quadMeshFile);
            Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> meshEntities = meshReader.ReadMesh();
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);
            foreach (var element in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, t);
                model.AddElement(new XContinuumElement2D(element.ElementType, element.Nodes, materialField,
                    integration, jIntegration));
            }
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(model, 1e-6);

            // Constraints
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            model.AddConstraint(bottomLeftNode, DisplacementDof.X, 0.0);
            if (fixBottom)
            {
                foreach (var node in finder.FindNodesWithY(0.0))
                {
                    model.AddConstraint(node, DisplacementDof.Y, 0.0);
                }
            }
            else
            {
                XNode2D bottomRightNode = finder.FindNodeWith(bottomWidth, 0.0);
                model.AddConstraint(bottomRightNode, DisplacementDof.Y, 0.0);
                model.AddConstraint(bottomLeftNode, DisplacementDof.Y, 0.0);
            }


            // Loads
            double distributedLoad = load / topWidth;
            var distrubutor = new LoadDistributor();
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(totalHeight);
            double[,] topLoads = distrubutor.DistributeLoad(topNodes, 0.0, distributedLoad);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                model.AddNodalLoad(topNodes[i], DisplacementDof.X, topLoads[i, 0]);
                model.AddNodalLoad(topNodes[i], DisplacementDof.Y, topLoads[i, 1]);
            }
        }

        private void HandleCrack()
        {
            var boundary = new FilletBoundary();
            var mesh = new BiMesh2D(model.Nodes, model.Elements, boundary);

            propagator = new Propagator(crack.Mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            crack = new TrackingExteriorCrackLSM(propagator);
            crack.Mesh = mesh;

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);

            // Mesh geometry interaction
            double webLeft = 0.5 * (bottomWidth - topWidth);
            var crackMouth = new CartesianPoint2D(webLeft, crackHeight);
            var crackTip = new CartesianPoint2D(webLeft + crackLength, crackHeight);
            crack.InitializeGeometry(crackMouth, crackTip);
        }

        private IReadOnlyList<ICartesianPoint2D> Analyze()
        {
            
            var solver = new SkylineSolver(model);
            QuasiStaticAnalysis analysis = new QuasiStaticAnalysis(model, crack.Mesh, crack, solver, fractureToughness, maxIterations);
             analysis.Analyze();
            return crack.CrackPath;
        }

        private class FilletBoundary: IDomainBoundary
        {
            private readonly double voidRectWidth, centerY, leftCenterX, rightCenterX;

            public FilletBoundary()
            {
                voidRectWidth = 0.5 * (bottomWidth - topWidth);
                centerY = flangeHeight + radius;
                leftCenterX = voidRectWidth - radius;
                rightCenterX = bottomWidth - voidRectWidth + radius;
            }

            public bool IsInside(ICartesianPoint2D point)
            {
                // Shapes
                var rectHull = new RectangularBoundary(0.0, bottomWidth, 0.0, totalHeight);
                var leftVoid = new RectangularBoundary(0.0, voidRectWidth, flangeHeight, totalHeight);
                var rightVoid = new RectangularBoundary(bottomWidth - voidRectWidth, bottomWidth, flangeHeight, totalHeight);
                var leftCircle = new Circle2D(new CartesianPoint2D(leftCenterX, centerY), radius);
                var rightCircle = new Circle2D(new CartesianPoint2D(rightCenterX, centerY), radius);

                if (rectHull.IsInside(point))
                {
                    if (leftVoid.IsInside(point)) // Over flange, left of web
                    {
                        if ((point.X > leftCenterX) && (point.Y < centerY))
                        {
                            if (leftCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside)
                            {
                                return true; // Inside left fillet
                            }
                            else return false;
                        }
                        else return false;
                    }
                    else if (rightVoid.IsInside(point)) // Over flange, right of web
                    {
                        if ((point.X < leftCenterX) && (point.Y < centerY))
                        {
                            if (leftCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside)
                            {
                                return true; // Inside right fillet
                            }
                            else return false;
                        }
                        else return false;
                    }
                    else return true; // Inside the flange or the web
                }
                else return false;
            }
        }
    }
}
