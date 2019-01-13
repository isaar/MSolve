using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;


namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class Fillet: IBenchmark
    {
        public static Builder SetupBenchmark(bool writePropagationPath, bool plotLSM)
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Meshes\fillet.msh";
            string propagationPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Propagation\crack_growth.txt";
            string plotPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Plots";
            string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Timing";
            //string meshPath = @"C:\Users\seraf\Desktop\GRACM\Fillet\Meshes\fillet.msh";
            //string propagationPath = @"C:\Users\seraf\Desktop\GRACM\Fillet\Propagation\crack_growth.txt";
            //string plotPath = @"C:\Users\seraf\Desktop\GRACM\Fillet\Plots";
            //string timingPath = @"C:\Users\seraf\Desktop\GRACM\Fillet\Timing";

            double growthLength = 5; // mm. Must be sufficiently larger than the element size.
            var builder = new Builder(meshPath, growthLength, timingPath, propagationPath);
            builder.WritePropagation = writePropagationPath;
            builder.HeavisideEnrichmentTolerance = 0.01;
            builder.RigidBCs = true;
            builder.NumSubdomains = 5;

            builder.LsmPlotDirectory = plotLSM ? plotPath: null;
            builder.MaxIterations = 13;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0; 

            return builder;
        }


        #region constants
        /// <summary>
        /// Thickness of the whole domain
        /// </summary>
        private const double t = 1.0; // mm

        /// <summary>
        /// Young's modulus
        /// </summary>
        private const double E = 2.1e12; // kN/mm^2

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3;

        /// <summary>
        /// The material used for the J-integral computation. It msut be stored separately from individual element materials.
        /// </summary>
        private static readonly HomogeneousElasticMaterial2D globalHomogeneousMaterial =
            HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

        /// <summary>
        /// The maximum value that the effective SIF can reach before collapse occurs.
        /// </summary>
        private const double fractureToughness = double.MaxValue;

        /// <summary>
        /// Magnitude in kN of "opening" load
        /// </summary>
        private const double load = 20.0; // kN


        // Geometry (m)
        private const double bottomWidth = 375.0, topWidth = 75.0, radius = 20.0;// mm
        private const double flangeHeight = 75.0, totalHeight = 150.0; // mm
        private const double crackHeight = flangeHeight + radius, crackLength = 5, webLeft = 0.5 * (bottomWidth - topWidth); //mm
        private const double infCrackHeight = 90.0, supCrackHeight = 105.0; //mm

        #endregion

        private readonly bool rigidBCs;

        private readonly double heavisideTol;

        /// <summary>
        /// The length by which the crack grows in each iteration.
        /// </summary>
        private readonly double growthLength;

        /// <summary>
        /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
        /// radius is larger than the length of the crack segments.
        /// </summary>
        private readonly double jIntegralRadiusOverElementSize;

        private readonly string meshPath;
        private readonly string lsmPlotDirectory;
        private readonly string propagationPath;
        private readonly bool writePropagation;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;
        private readonly int numSubdomains;

        private TrackingExteriorCrackLSM crack;
        private BiMesh2D mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private Fillet(string meshPath, double growthLength, double jIntegralRadiusOverElementSize,
             string lsmOutputDirectory, string propagationPath, bool writePropagation, int maxIterations,
              bool rigidBCs, double heavisideTol, int numSubdomains)
        {
            this.meshPath = meshPath;
            this.growthLength = growthLength;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.lsmPlotDirectory = lsmOutputDirectory;
            this.propagationPath = propagationPath;
            this.writePropagation = writePropagation;
            this.maxIterations = maxIterations;
            this.rigidBCs = rigidBCs;
            this.heavisideTol = heavisideTol;
            this.numSubdomains = numSubdomains;
        }

        /// <summary>
        /// The crack geometry description
        /// </summary>
        public ICrackDescription Crack { get { return crack; } }

        public IDomainDecomposer Decomposer { get; private set; }

        public IReadOnlyList<double> GrowthAngles { get; private set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public Model2D Model { get; private set; }

        public string Name { get { return "GRACM Fillet"; } }

        public string PlotDirectory { get { return lsmPlotDirectory; } }

        public Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> PossibleEnrichments { get; private set; }

        public void Analyze(ISolver solver)
        {
            var analysis = new QuasiStaticAnalysis(Model, mesh, crack, solver, fractureToughness, maxIterations);
            analysis.Analyze();

            // Write crack path
            Console.WriteLine("Crack path:");
            foreach (var point in crack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();

            // Write growth angles, lengths and SIFs if necessary
            if (writePropagation)
            {
                using (var writer = new StreamWriter(propagationPath))
                {
                    PropagationLogger logger = crack.CrackTipPropagators[crack.CrackTips[0]].Logger;
                    int numIterations = logger.GrowthAngles.Count;
                    writer.WriteLine(numIterations);
                    for (int i = 0; i < numIterations; ++i)
                    {
                        writer.Write(logger.GrowthAngles[i]);
                        writer.Write(" " + logger.GrowthLengths[i]);
                        writer.Write(" " + logger.SIFsMode1[i]);
                        writer.Write(" " + logger.SIFsMode2[i]);
                        writer.WriteLine();
                    }
                }
            }
        }

        public void InitializeModel()
        {
            Model = new Model2D();
            CreateMesh();
            ApplyBoundaryConditions();
            InitializeCrack();
            LimitEnrichedArea();
            DomainDecomposition();
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(Model, 1e-6);

            // Constraints
            if (rigidBCs)
            {
                foreach (var node in finder.FindNodesWithY(0.0))
                {
                    Model.AddConstraint(node, DisplacementDof.X, 0.0);
                    Model.AddConstraint(node, DisplacementDof.Y, 0.0);
                }
            }
            else // flexible
            {
                XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
                XNode2D bottomRightNode = finder.FindNodeWith(bottomWidth, 0.0);
                Model.AddConstraint(bottomLeftNode, DisplacementDof.X, 0.0);
                Model.AddConstraint(bottomLeftNode, DisplacementDof.Y, 0.0);
                Model.AddConstraint(bottomRightNode, DisplacementDof.X, 0.0);
                Model.AddConstraint(bottomRightNode, DisplacementDof.Y, 0.0);
            }

            // Loads
            double distributedLoad = load / topWidth;
            var distrubutor = new LoadDistributor();
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(totalHeight);
            double[,] topLoads = distrubutor.DistributeLoad(topNodes, 0.0, distributedLoad);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                Model.AddNodalLoad(topNodes[i], DisplacementDof.X, topLoads[i, 0]);
                Model.AddNodalLoad(topNodes[i], DisplacementDof.Y, topLoads[i, 1]);
            }
        }

        private void CreateMesh()
        {
            // Mesh generation
            var reader = new GmshReader(meshPath, true);
            Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> meshEntities = reader.ReadMesh();
            var nodes = new XNode2D[meshEntities.Item1.Count];
            var elementConnectivity = new List<XNode2D[]>();
            nodes = meshEntities.Item1.ToArray();
            foreach (var element in meshEntities.Item2)
            {
                elementConnectivity.Add(element.Nodes.ToArray());
            }

            // Nodes
            foreach (XNode2D node in nodes) Model.AddNode(node);

            // Integration rules
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            // Elements
            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                Model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }

            // Mesh usable for crack-mesh interaction
            
            mesh = new BiMesh2D(Model.Nodes, Model.Elements, new FilletBoundary());
        }

        private void DomainDecomposition() //TODO: this should not be hardcoded, but provided by the caller of the solver
        {
            //WARNING: Needs adjusting if the geometric constants change
            if (numSubdomains == 1)
            {
                var regions = new PolygonalRegion[1];
                var vertices = new CartesianPoint2D[4];
                vertices[0] = new CartesianPoint2D(0.0, 0.0);
                vertices[1] = new CartesianPoint2D(375.0, 0.0);
                vertices[2] = new CartesianPoint2D(375.0, 150.0);
                vertices[3] = new CartesianPoint2D(0.0, 150.0);
                var boundaries = new HashSet<LineSegment2D>();
                regions[0] = new PolygonalRegion(vertices, boundaries);
                Decomposer = new GuideDecomposer(regions, mesh);
            }
            else if (numSubdomains == 3)
            {
                double boundary1X = 176.0, boundary2X = 201.0;

                var regions = new PolygonalRegion[3];

                var vertices1 = new CartesianPoint2D[4];
                vertices1[0] = new CartesianPoint2D(0.0, 0.0);
                vertices1[1] = new CartesianPoint2D(boundary1X, 0.0);
                vertices1[2] = new CartesianPoint2D(boundary1X, 150.0);
                vertices1[3] = new CartesianPoint2D(0.0, 150.0);
                var boundaries1 = new HashSet<LineSegment2D>();
                boundaries1.Add(new LineSegment2D(vertices1[1], vertices1[2]));
                regions[0] = new PolygonalRegion(vertices1, boundaries1);

                var vertices2 = new CartesianPoint2D[4];
                vertices2[0] = new CartesianPoint2D(boundary1X, 0.0);
                vertices2[1] = new CartesianPoint2D(boundary2X, 0.0);
                vertices2[2] = new CartesianPoint2D(boundary2X, 150.0);
                vertices2[3] = new CartesianPoint2D(boundary1X, 150.0);
                var boundaries2 = new HashSet<LineSegment2D>();
                boundaries2.Add(new LineSegment2D(vertices2[3], vertices2[1]));
                boundaries2.Add(new LineSegment2D(vertices2[1], vertices2[2]));
                regions[1] = new PolygonalRegion(vertices2, boundaries2);

                var vertices3 = new CartesianPoint2D[4];
                vertices3[0] = new CartesianPoint2D(boundary2X, 0.0);
                vertices3[1] = new CartesianPoint2D(375.0, 0.0);
                vertices3[2] = new CartesianPoint2D(375.0, 150.0);
                vertices3[3] = new CartesianPoint2D(boundary2X, 150.0);
                var boundaries3 = new HashSet<LineSegment2D>();
                boundaries3.Add(new LineSegment2D(vertices3[3], vertices3[1]));
                regions[2] = new PolygonalRegion(vertices3, boundaries3);

                //Decomposer = new GuideDecomposer(regions, mesh);
                Decomposer = new TipAdaptiveDecomposer(mesh, regions, crack, new GuideDecomposer(regions, mesh));
            }
            else if (numSubdomains == 5)
            {
                double minX = 0, maxX = bottomWidth, minY = 0, maxY = totalHeight;
                double boundary1X = 165.0, boundary2X = 180.0, boundary3X = 195.0, boundary4X = 210.0;


                var regions = new PolygonalRegion[5];

                // Left crack
                regions[0] = DecomposeRectangle(minX, minY, boundary1X, maxY,
                    new bool[] { false, true, true, false });
                regions[1] = DecomposeRectangle(boundary1X, minY, boundary2X, maxY,
                    new bool[] { false, true, true, true });
                regions[2] = DecomposeRectangle(boundary2X, minY, boundary3X, maxY,
                    new bool[] { false, true, true, true });
                regions[3] = DecomposeRectangle(boundary3X, minY, boundary4X, maxY,
                    new bool[] { false, false, true, true });
                regions[4] = DecomposeRectangle(boundary4X, minY, maxX, maxY,
                    new bool[] { false, false, true, true });

                //Decomposer = new GuideDecomposer(regions, mesh);
                Decomposer = new TipAdaptiveDecomposer(mesh, regions, crack, new GuideDecomposer(regions, mesh));
            }
            else throw new NotImplementedException();
        }

        private PolygonalRegion DecomposeRectangle(double x1, double y1, double x2, double y2, bool[] areBoundaries)
        {
            var vertices = new CartesianPoint2D[4];
            vertices[0] = new CartesianPoint2D(x1, y1);
            vertices[1] = new CartesianPoint2D(x2, y1);
            vertices[2] = new CartesianPoint2D(x2, y2);
            vertices[3] = new CartesianPoint2D(x1, y2);
            var boundaries = new HashSet<LineSegment2D>();
            for (int i = 0; i < areBoundaries.Length; ++i)
            {
                if (areBoundaries[i])
                {
                    CartesianPoint2D start = vertices[i];
                    CartesianPoint2D end = vertices[(i + 1) % vertices.Length];
                    boundaries.Add(new LineSegment2D(start, end));
                }
            }
            return new PolygonalRegion(vertices, boundaries);
        }

        private void InitializeCrack()
        {
            IPropagator propagator;
            if (!writePropagation) propagator = new FixedPropagator(propagationPath, null);
            else
            {
                propagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            }

            var crackMouth = new CartesianPoint2D(webLeft, crackHeight);
            var crackTip = new CartesianPoint2D(webLeft + crackLength, crackHeight);
            var initialCrack = new PolyLine2D(crackMouth, crackTip);
            var lsmCrack = new TrackingExteriorCrackLSM(propagator, 0.0, new RelativeAreaResolver(heavisideTol));
            lsmCrack.Mesh = mesh;

            // Create enrichments          
            lsmCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(lsmCrack);
            lsmCrack.CrackTipEnrichments = new CrackTipEnrichments2D(lsmCrack, CrackTipPosition.Single);
            if (lsmPlotDirectory != null)
            {
                lsmCrack.EnrichmentLogger = new EnrichmentLogger(Model, lsmCrack, lsmPlotDirectory);
                lsmCrack.LevelSetLogger = new LevelSetLogger(Model, lsmCrack, lsmPlotDirectory);
                lsmCrack.LevelSetComparer = new PreviousLevelSetComparer(lsmCrack, lsmPlotDirectory);
            }

            // Mesh geometry interaction
            lsmCrack.InitializeGeometry(initialCrack);
            this.crack = lsmCrack;
        }

        private void LimitEnrichedArea()
        {
            var enrichedNodes = Model.Nodes.Where(node => (node.Y >= infCrackHeight) && (node.Y <= supCrackHeight)).ToList();
            PossibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();
            PossibleEnrichments.Add(crack.CrackBodyEnrichment, enrichedNodes);
            PossibleEnrichments.Add(crack.CrackTipEnrichments, enrichedNodes);
        }

        public class Builder: IBenchmarkBuilder
        {
            private readonly double growthLength;
            private readonly string meshPath;
            private readonly string propagationPath;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="meshPath">The absolute path of the mesh file.</param>
            /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
            /// <param name="timingDirectory">The absolute path of the file where slover timing will be written.</param>
            public Builder(string meshPath, double growthLength, string timingDirectory, string propagationPath)
            {
                this.growthLength = growthLength;
                this.meshPath = meshPath;
                this.TimingOutputDirectory = timingDirectory;
                this.propagationPath = propagationPath;
            }

            /// <summary>
            /// If set to true, the global y dof will be constrained for all nodes along the bottom edge, in order to simulate 
            /// a very thick rigid I-beam. If set to false, the global � degree of freedom is constrained only for the corner
            ///  nodes of the bottom edge, in order to simulate a very thin, flexible I-beam
            /// </summary>
            public bool RigidBCs { get; set; } = true;

            /// <summary>
            /// A node that lies in the positive halfplane defined by the body level set of a crack, will be enriched with 
            /// heaviside enrichment if Apos/(Apos+Aneg) &gt; <see cref="HeavisideEnrichmentTolerance"/> where Apos, Aneg 
            /// are the subareas of its nodal support  in the positive and negative halfplanes respectively. Similarly a
            /// node in the negative halfplane will be enriched if Aneg/(Apos+Aneg) &gt; 
            /// <see cref="HeavisideEnrichmentTolerance"/>.
            /// </summary>
            public double HeavisideEnrichmentTolerance { get; set; } = 0.0001;

            /// <summary>
            /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
            /// radius is larger than the length of the crack segments.
            /// </summary>
            public double JintegralRadiusOverElementSize { get; set; } = 2.0;

            /// <summary>
            /// The absolute path of the directory where output vtk files with the crack path and the level set functions at 
            /// each iteration will be written. Leave it null to avoid the performance cost it will introduce.
            /// </summary>
            public string LsmPlotDirectory { get; set; } = null;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = int.MaxValue;

            public int NumSubdomains { get; set; } = 1;

            /// <summary>
            /// The absolute path of the file where slover timing will be written.
            /// </summary>
            public string TimingOutputDirectory { get; }

            public bool WritePropagation { get; set; } = true;

            public IBenchmark BuildBenchmark()
            {
                return new Fillet(meshPath, growthLength, JintegralRadiusOverElementSize, LsmPlotDirectory,
                    propagationPath, WritePropagation, MaxIterations, RigidBCs, HeavisideEnrichmentTolerance, NumSubdomains);
            }
        }

        private class FilletBoundary : IDomainBoundary
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
