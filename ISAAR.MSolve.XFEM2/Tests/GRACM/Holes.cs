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
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Utilities;

//TODO: limit enriched nodes in the areas around the cracks
namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class Holes: IBenchmark
    {
        public static Builder SetupBenchmark(bool writePropagationPath, bool plotLSM)
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Meshes\holes.msh";
            string propagationPathLeft = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Propagation\left_growth.txt";
            string propagationPathRight = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Propagation\right_growth.txt";
            string plotPathLeft = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Plots\Left";
            string plotPathRight = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Plots\Right";
            string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Timing";
            //string meshPath = @"C:\Users\seraf\Desktop\GRACM\Holes\Meshes\holes.msh";
            //string propagationPathLeft = @"C:\Users\seraf\Desktop\GRACM\Holes\Propagation\left_growth.txt";
            //string propagationPathRight = @"C:\Users\seraf\Desktop\GRACM\Holes\Propagation\right_growth.txt";
            //string plotPathLeft = @"C:\Users\seraf\Desktop\GRACM\Holes\Plots\Left";
            //string plotPathRight = @"C:\Users\seraf\Desktop\GRACM\Holes\Plots\Right";
            //string timingPath = @"C:\Users\seraf\Desktop\GRACM\Holes\Timing";

            double growthLength = 1.0; // mm. Must be sufficiently larger than the element size. 0.8 works for left crack, 0.9 is so-so for right crack
            var builder = new Builder(meshPath, growthLength, timingPath, propagationPathLeft, propagationPathRight);
            builder.WritePropagation = writePropagationPath;

            builder.HeavisideEnrichmentTolerance = 0.12;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            // If you modify the following two parameters significantly, then you will need to redefine which nodes are expected 
            // to be enriched.
            builder.TipEnrichmentRadius = 0.0;
            builder.BC = BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY;
            builder.NumSubdomains = 8;

            builder.MaxIterations = 12;
            builder.LeftLsmPlotDirectory = plotLSM ? plotPathLeft : null;
            builder.RightLsmPlotDirectory = plotLSM ? plotPathRight : null;

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
        private const double E = 2e5; // N/mm^2

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3;

        /// <summary>
        /// Fracture toughness
        /// </summary>
        private const double KIc = double.MaxValue; // N.mm^(3/2) actually 1500

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
        /// Prescribed displacement
        /// </summary>
        private const double displacement = 0.1; //mm

        // Geometry
        private const double minX = 0.0, minY = 0.0, maxX = 20.0, maxY = 10.0; //mm
        private const double holeRadius = 2.0, leftHoleX = 3.0, leftHoleY = 7.0, rightHoleX = 17.0, rightHoleY = 3.0; //mm
        private const double initialCrackLength = 1.0;
        private const double leftCrackMouthX = 0.0, leftCrackMouthY = 2.85, leftCrackTipX = minX + initialCrackLength, leftCrackTipY = 2.85; //mm
        private const double rightCrackMouthX = 20.0, rightCrackMouthY = 7.15, rightCrackTipX = maxX - initialCrackLength, rightCrackTipY = 7.15; //mm
        #endregion

        private readonly BoundaryConditions bc;

        /// <summary>
        /// The length by which the crack grows in each iteration.
        /// </summary>
        private readonly double growthLength;

        private readonly double heavisideTol;

        /// <summary>
        /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
        /// radius is larger than the length of the crack segments.
        /// </summary>
        private readonly double jIntegralRadiusOverElementSize;

        private readonly string meshPath;
        private readonly string leftLsmPlotDirectory;
        private readonly string rightLsmPlotDirectory;
        private readonly string leftPropagationPath;
        private readonly string rightPropagationPath;
        private readonly bool writePropagation;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;
        private readonly int numSubdomains;
        private readonly double tipEnrichmentRadius;

        private TrackingExteriorCrackLSM leftCrack;
        private TrackingExteriorCrackLSM rightCrack;
        private BiMesh2D mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private Holes(string meshPath, double growthLength, BoundaryConditions bc, double jIntegralRadiusOverElementSize,
             double tipEnrichmentRadius, string leftLsmPlotDirectory, string rightLsmPlotDirectory,
             string leftPropagationPath, string rightPropagationPath, bool writePropagation, 
             int maxIterations, double heavisideTol, int numSubdomains)
        {
            this.meshPath = meshPath;
            this.growthLength = growthLength;
            this.bc = bc;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.tipEnrichmentRadius = tipEnrichmentRadius;
            this.leftLsmPlotDirectory = leftLsmPlotDirectory;
            this.rightLsmPlotDirectory = rightLsmPlotDirectory;
            this.leftPropagationPath = leftPropagationPath;
            this.rightPropagationPath = rightPropagationPath;
            this.writePropagation = writePropagation;
            this.maxIterations = maxIterations;
            this.heavisideTol = heavisideTol;
            this.numSubdomains = numSubdomains;
        }

        internal enum BoundaryConditions
        {
            BottomConstrainXY_TopDisplacementY, BottomConstrainXY_TopConstrainXDisplacementY,
            BottomConstrainY_TopDisplacementY, BottomDisplacementY_TopDisplacementY,
            BottomConstrainXDisplacementY_TopConstrainXDisplacementY
        }

        /// <summary>
        /// The crack geometry description
        /// </summary>
        public ICrackDescription Crack { get; private set; }

        public IDomainDecomposer Decomposer { get; private set; }

        public IReadOnlyList<XNode2D> EnrichedArea { get; private set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public Model2D Model { get; private set; }

        public string Name { get { return "GRACM Holes"; } }

        public string PlotDirectory { get { return leftLsmPlotDirectory; } }

        public Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> PossibleEnrichments { get; private set; }

        public void Analyze(ISolver solver)
        {
            var analysis = new QuasiStaticAnalysis(Model, mesh, Crack, solver, fractureToughness, maxIterations,
                new IntersectedMeshOutput(Model, Crack, leftLsmPlotDirectory + "\\field_output"));
            //var analysis = new QuasiStaticAnalysis(Model, mesh, Crack, solver, fractureToughness, maxIterations);
            analysis.Analyze();

            // Write crack path
            Console.WriteLine("Left crack path:");
            foreach (var point in leftCrack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();
            Console.WriteLine("Right crack path:");
            foreach (var point in rightCrack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();

            // Write growth angles, lengths and SIFs if necessary
            if (writePropagation)
            {
                using (var writer = new StreamWriter(leftPropagationPath))
                {
                    PropagationLogger logger = leftCrack.CrackTipPropagators[leftCrack.CrackTips[0]].Logger;
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
                using (var writer = new StreamWriter(rightPropagationPath))
                {
                    PropagationLogger logger = rightCrack.CrackTipPropagators[rightCrack.CrackTips[0]].Logger;
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
            XNode2D leftTopCorner = finder.FindNodeWith(minX, maxY);
            XNode2D rightBottomCorner = finder.FindNodeWith(maxX, minY);
            IReadOnlyList<XNode2D> bottomNodes = finder.FindNodesWithY(minY);
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(maxY);

            if (bc == BoundaryConditions.BottomConstrainXY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.X, 0.0);
                    Model.AddConstraint(node, DisplacementDof.Y, 0.0);
                }
                foreach (var node in topNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.Y, displacement);
                }
            }
            else if (bc == BoundaryConditions.BottomConstrainXY_TopConstrainXDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.X, 0.0);
                    Model.AddConstraint(node, DisplacementDof.Y, 0.0);
                }
                foreach (var node in topNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.X, 0.0);
                    Model.AddConstraint(node, DisplacementDof.Y, displacement);
                }
            }
            else if (bc == BoundaryConditions.BottomConstrainY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.Y, 0.0);
                }
                foreach (var node in topNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.Y, displacement);
                }             
                Model.AddConstraint(leftTopCorner, DisplacementDof.X, 0.0);
                Model.AddConstraint(rightBottomCorner, DisplacementDof.X, 0.0);
            }
            else if (bc == BoundaryConditions.BottomDisplacementY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.Y, - displacement / 2.0);
                }
                foreach (var node in topNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.Y, displacement / 2.0);
                }
                Model.AddConstraint(leftTopCorner, DisplacementDof.X, 0.0);
                Model.AddConstraint(rightBottomCorner, DisplacementDof.X, 0.0);
            }
            else if (bc == BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.X, 0.0);
                    Model.AddConstraint(node, DisplacementDof.Y, - displacement / 2.0);
                }
                foreach (var node in topNodes)
                {
                    Model.AddConstraint(node, DisplacementDof.X, 0.0);
                    Model.AddConstraint(node, DisplacementDof.Y, displacement / 2.0);
                }
            }
            else throw new Exception("Shouldn't have been reached.");
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

            mesh = new BiMesh2D(Model.Nodes, Model.Elements, new HolesBoundary());
        }

        private void DomainDecomposition() //TODO: this should not be hardcoded, but provided by the caller of the solver
        {
            //WARNING: Needs adjusting if the geometric constants change
            if (numSubdomains == 1)
            {
                var regions = new PolygonalRegion[1];
                regions[0] = DecomposeRectangle(minX, minY, maxX, maxY, new bool[] { false, false, false, false });
                Decomposer = new GuideDecomposer(regions, mesh);
            }
            else if (numSubdomains == 8)
            {
                var regions = new PolygonalRegion[8];
                double boundaryY = 5.0;
                double leftBoundary1X = 3.5, leftBoundary2X = 7.0, leftBoundary3X = 10.5;
                double rightBoundary1X = 9.5, rightBoundary2X = 13.0, rightBoundary3X = 16.5;

                // Left crack
                regions[0] = DecomposeRectangle(minX, minY, leftBoundary1X, boundaryY, 
                    new bool[] { false, true, true, false });
                regions[1] = DecomposeRectangle(leftBoundary1X, minY, leftBoundary2X, boundaryY, 
                    new bool[] { false, true, true, true });
                regions[2] = DecomposeRectangle(leftBoundary2X, minY, leftBoundary3X, boundaryY,
                    new bool[] { false, true, true, true });
                regions[3] = DecomposeRectangle(leftBoundary3X, minY, maxX, boundaryY,
                    new bool[] { false, false, true, true });

                // Right crack
                regions[4] = DecomposeRectangle(minX, boundaryY, rightBoundary1X, maxY,
                    new bool[] { true, true, false, false });
                regions[5] = DecomposeRectangle(rightBoundary1X, boundaryY, rightBoundary2X, maxY,
                    new bool[] { true, true, false, true });
                regions[6] = DecomposeRectangle(rightBoundary2X, boundaryY, rightBoundary3X, maxY,
                    new bool[] { true, true, false, true });
                regions[7] = DecomposeRectangle(rightBoundary3X, boundaryY, maxX, maxY,
                    new bool[] { true, false, false, true });

                //Decomposer = new GuideDecomposer(regions, mesh);
                Decomposer = new TipAdaptiveDecomposer(mesh, regions, Crack, new GuideDecomposer(regions, mesh));
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
            // Left crack
            // Propagator
            IPropagator leftPropagator;
            if (!writePropagation) leftPropagator = new FixedPropagator(leftPropagationPath, null);
            else
            {
                leftPropagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            }

            var initialLeftCrack = new PolyLine2D(new CartesianPoint2D(leftCrackMouthX, leftCrackMouthY),
                new CartesianPoint2D(leftCrackTipX, leftCrackTipY));
            leftCrack = new TrackingExteriorCrackLSM(leftPropagator, tipEnrichmentRadius, new RelativeAreaResolver(heavisideTol));
            leftCrack.Mesh = mesh;

            // Create enrichments          
            leftCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(leftCrack);
            leftCrack.CrackTipEnrichments = new CrackTipEnrichments2D(leftCrack, CrackTipPosition.Single);
            if (leftLsmPlotDirectory != null)
            {
                leftCrack.EnrichmentLogger = new EnrichmentLogger(Model, leftCrack, leftLsmPlotDirectory);
                leftCrack.LevelSetLogger = new LevelSetLogger(Model, leftCrack, leftLsmPlotDirectory);
                leftCrack.LevelSetComparer = new PreviousLevelSetComparer(leftCrack, leftLsmPlotDirectory);
            }

            // Mesh geometry interaction
            leftCrack.InitializeGeometry(initialLeftCrack);


            // Right crack
            // Propagator
            IPropagator rightPropagator;
            if (!writePropagation) rightPropagator = new FixedPropagator(rightPropagationPath, null);
            else
            {
                rightPropagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            }

            var initialRightCrack = new PolyLine2D(new CartesianPoint2D(rightCrackMouthX, rightCrackMouthY),
                new CartesianPoint2D(rightCrackTipX, rightCrackTipY));
            rightCrack = new TrackingExteriorCrackLSM(rightPropagator, tipEnrichmentRadius, new RelativeAreaResolver(heavisideTol));
            rightCrack.Mesh = mesh;

            // Create enrichments          
            rightCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(rightCrack);
            rightCrack.CrackTipEnrichments = new CrackTipEnrichments2D(rightCrack, CrackTipPosition.Single);
            if (rightLsmPlotDirectory != null)
            {
                rightCrack.EnrichmentLogger = new EnrichmentLogger(Model, rightCrack, rightLsmPlotDirectory);
                rightCrack.LevelSetLogger = new LevelSetLogger(Model, rightCrack, rightLsmPlotDirectory);
                rightCrack.LevelSetComparer = new PreviousLevelSetComparer(rightCrack, rightLsmPlotDirectory);
            }

            // Mesh geometry interaction
            rightCrack.InitializeGeometry(initialRightCrack);


            // Container for both cracks
            Crack = new MultipleCracksDisjoint(new TrackingExteriorCrackLSM[] { leftCrack, rightCrack });
            //Crack = new MultipleCracksDisjoint(new TrackingExteriorCrackLSM[] { leftCrack });
        }

        private void LimitEnrichedArea()
        {
            // TODO: these need to be adjucted if the dimensions of the benchmark change.
            double leftMinX = 0, leftMaxX = 15.6, leftMinY = 1.6, leftMaxY = 4.2;
            double rightMinX = 4.5, rightMaxX = 20.0, rightMinY = 5.5, rightMaxY = 8.5;

            //TODO: do that
            //EnrichedArea = Model.Nodes.Where(node => (node.Y >= infCrackHeight) && (node.Y <= supCrackHeight)).ToList();
            PossibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();

            //var leftEnrichedNodes = Model.Nodes;
            var leftEnrichedNodes = Model.Nodes.Where(node =>
                (node.X >= leftMinX) && (node.X <= leftMaxX) && (node.Y >= leftMinY) && (node.Y <= leftMaxY)).ToList();
            PossibleEnrichments.Add(leftCrack.CrackBodyEnrichment, leftEnrichedNodes);
            PossibleEnrichments.Add(leftCrack.CrackTipEnrichments, leftEnrichedNodes);

            //var rightEnrichedNodes = Model.Nodes;
            var rightEnrichedNodes = Model.Nodes.Where(node =>
                (node.X >= rightMinX) && (node.X <= rightMaxX) && (node.Y >= rightMinY) && (node.Y <= rightMaxY)).ToList();
            PossibleEnrichments.Add(rightCrack.CrackBodyEnrichment, rightEnrichedNodes);
            PossibleEnrichments.Add(rightCrack.CrackTipEnrichments, rightEnrichedNodes);
        }

        public class Builder : IBenchmarkBuilder
        {
            private readonly double growthLength;
            private readonly string meshPath;
            private readonly string leftPropagationPath;
            private readonly string rightPropagationPath;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="meshPath">The absolute path of the mesh file.</param>
            /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
            /// <param name="timingDirectory">The absolute path of the file where slover timing will be written.</param>
            public Builder(string meshPath, double growthLength, string timingDirectory, 
                string leftPropagationPath, string rightPropagationPath)
            {
                this.growthLength = growthLength;
                this.meshPath = meshPath;
                this.TimingOutputDirectory = timingDirectory;
                this.leftPropagationPath = leftPropagationPath;
                this.rightPropagationPath = rightPropagationPath;
            }

            public BoundaryConditions BC { get; set; } = BoundaryConditions.BottomConstrainXY_TopDisplacementY;

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
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            public PropagationLogger KnownLeftPropagation { get; set; } = null;

            /// <summary>
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            public PropagationLogger KnownRightPropagation { get; set; } = null;

            /// <summary>
            /// The absolute path of the directory where output vtk files with the crack path and the level set functions at 
            /// each iteration will be written. Leave it null to avoid the performance cost it will introduce.
            /// </summary>
            public string LeftLsmPlotDirectory { get; set; } = null;

            /// <summary>
            /// The absolute path of the directory where output vtk files with the crack path and the level set functions at 
            /// each iteration will be written. Leave it null to avoid the performance cost it will introduce.
            /// </summary>
            public string RightLsmPlotDirectory { get; set; } = null;

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

            public double TipEnrichmentRadius { get; set; } = 0.0;

            public bool WritePropagation { get; set; } = true;

            public IBenchmark BuildBenchmark()
            {
                return new Holes(meshPath, growthLength, BC, JintegralRadiusOverElementSize, TipEnrichmentRadius,
                    LeftLsmPlotDirectory, RightLsmPlotDirectory, leftPropagationPath, rightPropagationPath, WritePropagation, 
                    MaxIterations, HeavisideEnrichmentTolerance, NumSubdomains);
            }
        }

        private class HolesBoundary : IDomainBoundary
        {
            public HolesBoundary()
            {
            }

            public bool IsInside(ICartesianPoint2D point)
            {
                // Shapes
                var rectHull = new RectangularBoundary(minX, maxX, minY, maxY);
                var leftCircle = new Circle2D(new CartesianPoint2D(leftHoleX, leftHoleY), holeRadius);
                var rightCircle = new Circle2D(new CartesianPoint2D(rightHoleX, rightHoleY), holeRadius);

                // Intrnal points lie inside the rectangle, but outside the circular holes.
                if (rectHull.IsInside(point))
                {
                    if (leftCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside) return true;
                    if (rightCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside) return true;
                }
                return false;
            }
        }
    }
}
