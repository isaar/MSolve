using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.GMSH;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.XFEM.Analyzers;
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
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests
{
    public class HolesBenchmark
    {
        public enum BoundaryConditions
        {
            BottomConstrainXY_TopDisplacementY, BottomConstrainXY_TopConstrainXDisplacementY,
            BottomConstrainY_TopDisplacementY, BottomDisplacementY_TopDisplacementY,
            BottomConstrainXDisplacementY_TopConstrainXDisplacementY
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
        
        private const int subdomainID = 0;
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
        private readonly string subdomainPlotDirectory;
        private readonly string leftPropagationPath;
        private readonly string rightPropagationPath;
        private readonly bool writePropagation;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;
        private readonly double tipEnrichmentRadius;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private HolesBenchmark(string meshPath, double growthLength, BoundaryConditions bc, double jIntegralRadiusOverElementSize,
             double tipEnrichmentRadius, string leftLsmPlotDirectory, string rightLsmPlotDirectory, string subdomainPlotDirectory,
             string leftPropagationPath, string rightPropagationPath, bool writePropagation,
             int maxIterations, double heavisideTol)
        {
            this.meshPath = meshPath;
            this.growthLength = growthLength;
            this.bc = bc;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.tipEnrichmentRadius = tipEnrichmentRadius;
            this.leftLsmPlotDirectory = leftLsmPlotDirectory;
            this.rightLsmPlotDirectory = rightLsmPlotDirectory;
            this.subdomainPlotDirectory = subdomainPlotDirectory;
            this.leftPropagationPath = leftPropagationPath;
            this.rightPropagationPath = rightPropagationPath;
            this.writePropagation = writePropagation;
            this.maxIterations = maxIterations;
            this.heavisideTol = heavisideTol;
        }

        /// <summary>
        /// The crack geometry description
        /// </summary>
        public ICrackDescription Crack { get; private set; }

        public BidirectionalMesh2D<XNode, XContinuumElement2D> Mesh { get; set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public XModel Model { get; private set; }

        public string Name { get { return "Twin Holes benchmark"; } }

        public string PlotDirectory { get { return leftLsmPlotDirectory; } }

        public TrackingExteriorCrackLSM LeftCrack { get; set; }
        public TrackingExteriorCrackLSM RightCrack { get; set; }

        public void Analyze(ISolver solver)
        {
            var problem = new ProblemStructural(Model, solver);
            var analyzer = new QuasiStaticCrackPropagationAnalyzer(Model, solver, /*problem,*/ Crack, fractureToughness,
                maxIterations);

            // Subdomain plots
            if (subdomainPlotDirectory != null)
            {
                if (solver is FetiDPSolver fetiDP)
                {
                    analyzer.DDLogger = new DomainDecompositionLoggerFetiDP(fetiDP, subdomainPlotDirectory);
                }
                else analyzer.DDLogger = new DomainDecompositionLogger(subdomainPlotDirectory);
            }

            analyzer.Initialize();
            analyzer.Analyze();

            #region crack propagation output
            //// Write crack path
            //Console.WriteLine("Left crack path:");
            //foreach (var point in leftCrack.CrackPath)
            //{
            //    Console.WriteLine("{0} {1}", point.X, point.Y);
            //}
            //Console.WriteLine();
            //Console.WriteLine("Right crack path:");
            //foreach (var point in rightCrack.CrackPath)
            //{
            //    Console.WriteLine("{0} {1}", point.X, point.Y);
            //}
            //Console.WriteLine();

            //// Write growth angles, lengths and SIFs if necessary
            //if (writePropagation)
            //{
            //    using (var writer = new StreamWriter(leftPropagationPath))
            //    {
            //        PropagationLogger logger = leftCrack.CrackTipPropagators[leftCrack.CrackTips[0]].Logger;
            //        int numIterations = logger.GrowthAngles.Count;
            //        writer.WriteLine(numIterations);
            //        for (int i = 0; i < numIterations; ++i)
            //        {
            //            writer.Write(logger.GrowthAngles[i]);
            //            writer.Write(" " + logger.GrowthLengths[i]);
            //            writer.Write(" " + logger.SIFsMode1[i]);
            //            writer.Write(" " + logger.SIFsMode2[i]);
            //            writer.WriteLine();
            //        }
            //    }
            //    using (var writer = new StreamWriter(rightPropagationPath))
            //    {
            //        PropagationLogger logger = rightCrack.CrackTipPropagators[rightCrack.CrackTips[0]].Logger;
            //        int numIterations = logger.GrowthAngles.Count;
            //        writer.WriteLine(numIterations);
            //        for (int i = 0; i < numIterations; ++i)
            //        {
            //            writer.Write(logger.GrowthAngles[i]);
            //            writer.Write(" " + logger.GrowthLengths[i]);
            //            writer.Write(" " + logger.SIFsMode1[i]);
            //            writer.Write(" " + logger.SIFsMode2[i]);
            //            writer.WriteLine();
            //        }
            //    }
            //}
            #endregion
        }

        public void InitializeModel()
        {
            Model = new XModel();
            CreateMesh();
            ApplyBoundaryConditions();
            InitializeCrack();
        }

        private void ApplyBoundaryConditions()
        {
            double meshTol = 1E-6;
            XNode leftTopCorner = Model.Nodes.Where(
                node => (Math.Abs(node.X - minX) <= meshTol) && (Math.Abs(node.Y - maxY) <= meshTol)).First();
            XNode rightBottomCorner = Model.Nodes.Where(
                node => (Math.Abs(node.X - maxX) <= meshTol) && (Math.Abs(node.Y - minY) <= meshTol)).First();
            XNode[] bottomNodes = Model.Nodes.Where(node => Math.Abs(node.Y - minY) <= meshTol).ToArray();
            XNode[] topNodes = Model.Nodes.Where(node => Math.Abs(node.Y - maxY) <= meshTol).ToArray();

            if (bc == BoundaryConditions.BottomConstrainXY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = displacement });
                }
            }
            else if (bc == BoundaryConditions.BottomConstrainXY_TopConstrainXDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = displacement });
                }
            }
            else if (bc == BoundaryConditions.BottomConstrainY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = displacement });
                }
                leftTopCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                rightBottomCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            }
            else if (bc == BoundaryConditions.BottomDisplacementY_TopDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = -0.5 * displacement });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.5 * displacement });
                }
                leftTopCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                rightBottomCorner.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            }
            else if (bc == BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY)
            {
                foreach (var node in bottomNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = -0.5 * displacement });
                }
                foreach (var node in topNodes)
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.5 * displacement });
                }
            }
            else throw new Exception("Shouldn't have been reached.");
        }

        private void CreateMesh()
        {
            Model.Subdomains.Add(subdomainID, new XSubdomain(subdomainID));

            // Mesh generation
            var reader = new GmshReader<XNode>(meshPath);
            (IReadOnlyList<XNode> nodes, IReadOnlyList<CellConnectivity<XNode>> elementConnectivities) = reader.CreateMesh(
                (id, x, y, z) => new XNode(id, x, y, z));

            // Nodes
            foreach (XNode node in nodes) Model.Nodes.Add(node);

            // Integration rules
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)));
            var jIntegration =
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(4, 4));

            // Elements
            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
            var factory = new XContinuumElement2DFactory(integration, jIntegration, material);
            var cells = new XContinuumElement2D[elementConnectivities.Count];
            for (int e = 0; e < cells.Length; ++e)
            {
                XContinuumElement2D element = factory.CreateElement(e, CellType.Quad4, elementConnectivities[e].Vertices);
                cells[e] = element;
                Model.Elements.Add(element);
                Model.Subdomains[subdomainID].Elements.Add(Model.Elements[e]);
            }

            // Mesh usable for crack-mesh interaction
            var boundary = new HolesBoundary();
            Model.Boundary = boundary;
            Mesh = new BidirectionalMesh2D<XNode, XContinuumElement2D>(Model.Nodes, cells, boundary);
        }

        private void InitializeCrack()
        {
            // Left crack
            // Propagator
            IPropagator leftPropagator;
            if (!writePropagation) leftPropagator = new FixedPropagator(leftPropagationPath, null);
            else
            {
                leftPropagator = new Propagator(Mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            }

            var initialLeftCrack = new PolyLine2D(new CartesianPoint(leftCrackMouthX, leftCrackMouthY),
                new CartesianPoint(leftCrackTipX, leftCrackTipY));
            LeftCrack = new TrackingExteriorCrackLSM(leftPropagator, tipEnrichmentRadius, new RelativeAreaResolver(heavisideTol));
            LeftCrack.Mesh = Mesh;

            // Create enrichments          
            LeftCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(LeftCrack);
            LeftCrack.CrackTipEnrichments = new CrackTipEnrichments2D(LeftCrack, CrackTipPosition.Single);
            if (leftLsmPlotDirectory != null)
            {
                LeftCrack.EnrichmentLogger = new EnrichmentLogger(Model, LeftCrack, leftLsmPlotDirectory);
                LeftCrack.LevelSetLogger = new LevelSetLogger(Model, LeftCrack, leftLsmPlotDirectory);
                LeftCrack.LevelSetComparer = new PreviousLevelSetComparer(LeftCrack, leftLsmPlotDirectory);
            }

            // Mesh geometry interaction
            LeftCrack.InitializeGeometry(initialLeftCrack);

            // Right crack
            // Propagator
            IPropagator rightPropagator;
            if (!writePropagation) rightPropagator = new FixedPropagator(rightPropagationPath, null);
            else
            {
                rightPropagator = new Propagator(Mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            }

            var initialRightCrack = new PolyLine2D(new CartesianPoint(rightCrackMouthX, rightCrackMouthY),
                new CartesianPoint(rightCrackTipX, rightCrackTipY));
            RightCrack = new TrackingExteriorCrackLSM(rightPropagator, tipEnrichmentRadius, new RelativeAreaResolver(heavisideTol));
            RightCrack.Mesh = Mesh;

            // Create enrichments          
            RightCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(RightCrack);
            RightCrack.CrackTipEnrichments = new CrackTipEnrichments2D(RightCrack, CrackTipPosition.Single);
            if (rightLsmPlotDirectory != null)
            {
                RightCrack.EnrichmentLogger = new EnrichmentLogger(Model, RightCrack, rightLsmPlotDirectory);
                RightCrack.LevelSetLogger = new LevelSetLogger(Model, RightCrack, rightLsmPlotDirectory);
                RightCrack.LevelSetComparer = new PreviousLevelSetComparer(RightCrack, rightLsmPlotDirectory);
            }

            // Mesh geometry interaction
            RightCrack.InitializeGeometry(initialRightCrack);

            // Container for both cracks
            Crack = new MultipleCracksDisjoint(new TrackingExteriorCrackLSM[] { LeftCrack, RightCrack });
            //Crack = new MultipleCracksDisjoint(new TrackingExteriorCrackLSM[] { leftCrack });
        }

        public class Builder //: IBenchmarkBuilder
        {
            private readonly double growthLength;
            private readonly string meshPath;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="meshPath">The absolute path of the mesh file.</param>
            /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
            /// <param name="timingDirectory">The absolute path of the file where slover timing will be written.</param>
            public Builder(string meshPath, double growthLength)
            {
                this.growthLength = growthLength;
                this.meshPath = meshPath;
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

            public string LeftPropagationPath { get; set; }
            public string RightPropagationPath { get; set; }

            /// <summary>
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            //public PropagationLogger KnownLeftPropagation { get; set; } = null;

            /// <summary>
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            //public PropagationLogger KnownRightPropagation { get; set; } = null;

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

            public string SubdomainPlotDirectory { get; set; } = null;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = int.MaxValue;

            //public int NumSubdomains { get; set; } = 1;

            /// <summary>
            /// The absolute path of the file where slover timing will be written.
            /// </summary>
            public string TimingOutputDirectory { get; set; }

            public double TipEnrichmentRadius { get; set; } = 0.0;

            public bool WritePropagation { get; set; } = true;

            public HolesBenchmark BuildBenchmark()
            {
                return new HolesBenchmark(meshPath, growthLength, BC, JintegralRadiusOverElementSize, TipEnrichmentRadius,
                    LeftLsmPlotDirectory, RightLsmPlotDirectory, SubdomainPlotDirectory,
                    LeftPropagationPath, RightPropagationPath, WritePropagation,
                    MaxIterations, HeavisideEnrichmentTolerance);
            }
        }

        private class HolesBoundary : IDomain2DBoundary
        {
            public HolesBoundary()
            {
            }

            public bool IsInside(CartesianPoint point)
            {
                // Shapes
                var rectHull = new Rectangular2DBoundary(minX, maxX, minY, maxY);
                var leftCircle = new Circle2D(new CartesianPoint(leftHoleX, leftHoleY), holeRadius);
                var rightCircle = new Circle2D(new CartesianPoint(rightHoleX, rightHoleY), holeRadius);

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