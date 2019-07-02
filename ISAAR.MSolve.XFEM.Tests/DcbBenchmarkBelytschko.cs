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
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition;
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
    public class DcbBenchmarkBelytschko //: IBenchmark
    {
        #region constants
        ///// <summary>
        ///// The material used for the J-integral computation. It msut be stored separately from individual element materials.
        ///// </summary>
        //private static readonly HomogeneousElasticMaterial2D globalHomogeneousMaterial =
        //    HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

        /// <summary>
        /// The maximum value that the effective SIF can reach before collapse occurs.
        /// </summary>
        private const double fractureToughness = double.MaxValue;

        public static readonly double h = 3.94, L = 3 * h;// in
        private static readonly double v = 0.3, E = 3e7; // psi=lbs/in^2
        private static readonly double load = 197; // lbs
        private static readonly double a = 3.95, da = 0.5; // in 
        private static readonly double dTheta = 5.71 * Math.PI / 180; // initial crack angle
        #endregion

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

        private readonly int numElementsY;
        private readonly int numSubdomainsX;
        private readonly int numSubdomainsY;

        private readonly string lsmPlotDirectory;
        private readonly string subdomainPlotDirectory;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;

        private readonly double tipEnrichmentRadius = 0.0;

        private TrackingExteriorCrackLSM crack;
        private BidirectionalMesh2D<XNode, XContinuumElement2D> mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        public DcbBenchmarkBelytschko(int numElementsY, int numSubdomainsX, int numSubdomainsY, double growthLength, 
            double tipEnrichmentRadius, double jIntegralRadiusOverElementSize,  int maxIterations, double heavisideTol,
            string lsmPlotDirectory, string subdomainPlotDirectory)
        {
            this.numElementsY = numElementsY;
            this.numSubdomainsX = numSubdomainsX;
            this.numSubdomainsY = numSubdomainsY;
            this.growthLength = growthLength;
            this.tipEnrichmentRadius = tipEnrichmentRadius;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.lsmPlotDirectory = lsmPlotDirectory;
            this.subdomainPlotDirectory = subdomainPlotDirectory;
            this.maxIterations = maxIterations;
            this.heavisideTol = heavisideTol;
        }

        /// <summary>
        /// The crack geometry description. Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public TrackingExteriorCrackLSM Crack { get { return crack; } }

        public CartesianPoint CrackMouth { get; private set; }

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
                    analyzer.DDLogger = new DomainDecompositionLoggerFetiDP(fetiDP, subdomainPlotDirectory, true);
                }
                else analyzer.DDLogger = new DomainDecompositionLogger(subdomainPlotDirectory);
            }

            analyzer.Initialize();
            analyzer.Analyze();
        }

        public void InitializeModel()
        {
            Model = new XModel();
            CreateModel();
            InitializeCrack();
        }

        
        private void CreateModel()
        {
            var builder = new Uniform2DXModelBuilder();
            builder.DomainLengthX = L;
            builder.DomainLengthY = h;
            builder.NumSubdomainsX = numSubdomainsX;
            builder.NumSubdomainsY = numSubdomainsY;
            builder.NumTotalElementsX = 3 * numElementsY;
            builder.NumTotalElementsY = numElementsY;
            builder.YoungModulus = E;
            builder.PrescribeDisplacement(Uniform2DXModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DXModelBuilder.BoundaryRegion.RightSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DXModelBuilder.BoundaryRegion.UpperLeftCorner, StructuralDof.TranslationY, load);
            builder.DistributeLoadAtNodes(Uniform2DXModelBuilder.BoundaryRegion.LowerLeftCorner, StructuralDof.TranslationY, -load);

            (Model, mesh) = builder.BuildModel();
        }

        private void InitializeCrack()
        {
            var globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
            IPropagator propagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            CrackMouth = new CartesianPoint(0.0, h/2);
            var crackKink = new CartesianPoint(a, h / 2);
            var initialCrack = new PolyLine2D(CrackMouth, crackKink);
            initialCrack.UpdateGeometry(-dTheta, da);
            //var crackTip = new CartesianPoint(a + da * Math.Cos(dTheta), h/2 - da * Math.Sin(dTheta));

            var lsmCrack = new TrackingExteriorCrackLSM(propagator, tipEnrichmentRadius, new RelativeAreaResolver(heavisideTol));
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
            //lsmCrack.UpdateGeometry(-dTheta, da);
            this.crack = lsmCrack;
        }

        public class Builder //: IBenchmarkBuilder
        {
            private readonly int numElementsY;
            private readonly int numSubdomainsX;
            private readonly int numSubdomainsY;

            public Builder(int numElementsY, int numSubdomainsX, int numSubdomainsY)
            {
                this.numElementsY = numElementsY;
                this.numSubdomainsX = numSubdomainsX;
                this.numSubdomainsY = numSubdomainsY;
            }

            /// <summary>
            /// A node that lies in the positive halfplane defined by the body level set of a crack, will be enriched with 
            /// heaviside enrichment if Apos/(Apos+Aneg) &gt; <see cref="HeavisideEnrichmentTolerance"/> where Apos, Aneg 
            /// are the subareas of its nodal support  in the positive and negative halfplanes respectively. Similarly a
            /// node in the negative halfplane will be enriched if Aneg/(Apos+Aneg) &gt; 
            /// <see cref="HeavisideEnrichmentTolerance"/>.
            /// </summary>
            public double HeavisideEnrichmentTolerance { get; set; } = 0.0001;

            /// <summary>
            /// Must be sufficiently larger than the element size.
            /// </summary>
            public double GrowthLength { get; set; } = 0.3; //in

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
            public int MaxIterations { get; set; } = 8; //TODO: After that I noticed very weird behaviour

            public string PropagationPath { get; set; }

            /// <summary>
            /// The absolute path of the file where solver timing will be written.
            /// </summary>
            public string TimingOutputDirectory { get; }

            public double TipEnrichmentRadius { get; set; } = 0.0;

            public DcbBenchmarkBelytschko BuildBenchmark()
            {
                return new DcbBenchmarkBelytschko(numElementsY, numSubdomainsX, numSubdomainsY, GrowthLength, TipEnrichmentRadius, 
                    JintegralRadiusOverElementSize, MaxIterations, HeavisideEnrichmentTolerance,
                    LsmPlotDirectory, SubdomainPlotDirectory);
            }
        }
    }
}
