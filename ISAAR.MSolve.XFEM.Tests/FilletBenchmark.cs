using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.XFEM.Analyzers;
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
    public class FilletBenchmark //: IBenchmark
    {
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

        private const int subdomainID = 0;
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
        private readonly string subdomainPlotDirectory;

        private readonly string propagationPath;
        private readonly bool writePropagation;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;

        private TrackingExteriorCrackLSM crack;
        private BidirectionalMesh2D<XNode, XContinuumElement2D> mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        public FilletBenchmark(double growthLength, double jIntegralRadiusOverElementSize, string meshPath,
             string lsmPlotDirectory, string subdomainPlotDirectory, string propagationPath, bool writePropagation,
             int maxIterations, bool rigidBCs, double heavisideTol)
        {
            this.growthLength = growthLength;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.meshPath = meshPath;
            this.lsmPlotDirectory = lsmPlotDirectory;
            this.subdomainPlotDirectory = subdomainPlotDirectory;
            this.propagationPath = propagationPath;
            this.writePropagation = writePropagation;
            this.maxIterations = maxIterations;
            this.rigidBCs = rigidBCs;
            this.heavisideTol = heavisideTol;
        }

        /// <summary>
        /// The crack geometry description. Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public TrackingExteriorCrackLSM Crack { get { return crack; } }

        //public IReadOnlyList<double> GrowthAngles { get; private set; }

        public IMesh2D<XNode, XContinuumElement2D> Mesh => mesh;

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public XModel Model { get; private set; }

        public string Name { get { return "Fillet benchmark"; } }

        //public string PlotDirectory { get { return lsmPlotDirectory; } }

        public void Analyze(ISolver solver)
        {
            var problem = new ProblemStructural(Model, solver);
            var analyzer = new QuasiStaticCrackPropagationAnalyzer(Model, solver, /*problem,*/ crack, fractureToughness, 
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
            // Write crack path
            //Console.WriteLine("Crack path:");
            //foreach (var point in crack.CrackPath)
            //{
            //    Console.WriteLine("{0} {1}", point.X, point.Y);
            //}
            //Console.WriteLine();

            // Write growth angles, lengths and SIFs if necessary
            //if (writePropagation)
            //{
            //    using (var writer = new StreamWriter(propagationPath))
            //    {
            //        PropagationLogger logger = crack.CrackTipPropagators[crack.CrackTips[0]].Logger;
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
            CreateModel();
            ApplyBoundaryConditions();
            InitializeCrack();
        }

        private void ApplyBoundaryConditions()
        {
            double meshTol = 1E-6;
            // Constraints
            if (rigidBCs)
            {
                foreach (var node in Model.Nodes.Where(n => Math.Abs(n.Y) <= meshTol))
                {
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                    node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                }
            }
            else // flexible
            {
                XNode bottomLeftNode = Model.Nodes.Where(n => (Math.Abs(n.X) <= meshTol) && (Math.Abs(n.Y) <= meshTol)).First();
                XNode bottomRightNode = Model.Nodes.Where(
                    n => (Math.Abs(n.X - bottomWidth) <= meshTol) && (Math.Abs(n.Y) <= meshTol)).First();
                bottomLeftNode.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                bottomLeftNode.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
                bottomRightNode.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                bottomRightNode.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
            }

            // Distribute load amongst top nodes uniformly
            XNode[] topNodes = Model.Nodes.Where(n => Math.Abs(n.Y - totalHeight) <= meshTol).ToArray();
            //double distributedLoad = load / topWidth;
            for (int i = 0; i < topNodes.Length; ++i)
            {
                Model.NodalLoads.Add(new NodalLoad(topNodes[i], StructuralDof.TranslationY, load / topNodes.Length));
            }
        }

        private void CreateModel()
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
            var boundary = new FilletBoundary();
            Model.Boundary = boundary;
            mesh = new BidirectionalMesh2D<XNode, XContinuumElement2D>(Model.Nodes, cells, boundary);
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

            var crackMouth = new CartesianPoint(webLeft, crackHeight);
            var crackTip = new CartesianPoint(webLeft + crackLength, crackHeight);
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
            public Builder(string meshPath, double growthLength /*, string timingDirectory, string propagationPath*/)
            {
                this.meshPath = meshPath;
                this.growthLength = growthLength;
                //this.TimingOutputDirectory = timingDirectory;
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

            public string SubdomainPlotDirectory { get; set; } = null;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = int.MaxValue;

            //public int NumSubdomains { get; set; } = 1;

            public string PropagationPath { get; set; }

            /// <summary>
            /// The absolute path of the file where solver timing will be written.
            /// </summary>
            public string TimingOutputDirectory { get; }

            public bool WritePropagation { get; set; } = true;

            public FilletBenchmark BuildBenchmark()
            {
                return new FilletBenchmark(growthLength, JintegralRadiusOverElementSize, meshPath, LsmPlotDirectory,
                    SubdomainPlotDirectory, PropagationPath, WritePropagation, 
                    MaxIterations, RigidBCs, HeavisideEnrichmentTolerance);
            }
        }

        private class FilletBoundary : IDomain2DBoundary
        {
            private readonly double voidRectWidth, centerY, leftCenterX, rightCenterX;

            public FilletBoundary()
            {
                voidRectWidth = 0.5 * (bottomWidth - topWidth);
                centerY = flangeHeight + radius;
                leftCenterX = voidRectWidth - radius;
                rightCenterX = bottomWidth - voidRectWidth + radius;
            }

            public bool IsInside(CartesianPoint point)
            {
                // Shapes
                var rectHull = new Rectangular2DBoundary(0.0, bottomWidth, 0.0, totalHeight);
                var leftVoid = new Rectangular2DBoundary(0.0, voidRectWidth, flangeHeight, totalHeight);
                var rightVoid = new Rectangular2DBoundary(bottomWidth - voidRectWidth, bottomWidth, flangeHeight, totalHeight);
                var leftCircle = new Circle2D(new CartesianPoint(leftCenterX, centerY), radius);
                var rightCircle = new Circle2D(new CartesianPoint(rightCenterX, centerY), radius);

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
