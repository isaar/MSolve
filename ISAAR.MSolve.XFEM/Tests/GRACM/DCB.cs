using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class DCB
    {
        #region dimensions and material constants
        /// <summary>
        /// Horizontal dimension of the whole domain in inches.
        /// </summary>
        private const double L = 11.82;

        /// <summary>
        /// Vertical dimension of the whole domain in inches.
        /// </summary>
        private const double h = 3.94;

        /// <summary>
        /// Length of first crack segment (horizontal) in inches.
        /// </summary>
        private const double a = 3.95;

        /// <summary>
        /// Length of second crack segment (inclined) in inches.
        /// </summary>
        private const double da = 0.5;

        /// <summary>
        /// Initial crack angle in radians.
        /// </summary>
        private const double dTheta = 5.71 * Math.PI / 180;

        /// <summary>
        /// Thickness of the whole domain in inches.
        /// </summary>
        private const double t = 1.0;

        /// <summary>
        /// Young's modulus in psi=lbs/in^2
        /// </summary>
        private const double E = 3e7;

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3; // psi=lbs/in^2

        /// <summary>
        /// The material used for the J-integral computation. It msut be stored separately from individual element materials.
        /// </summary>
        private static readonly HomogeneousElasticMaterial2D globalHomogeneousMaterial = 
            HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);

        /// <summary>
        /// The maximum value that the effective SIF can reach before collapse occurs.
        /// </summary>
        private const double fractureToughness = double.MaxValue;

        /// <summary>
        /// Magnitude in lbs of "opening" point loads applied at the upper left and lower left corners of the domain.
        /// </summary>
        private const double load = 197; // lbs
        #endregion

        /// <summary>
        /// The approximate size of the elements. If the mesh is not uniform, this applies to the finest region.
        /// </summary>
        private readonly double elementSize;

        /// <summary>
        /// The length by which the crack grows in each iteration.
        /// </summary>
        private readonly double growthLength;

        /// <summary>
        /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
        /// radius is larger than the length of the crack segments.
        /// </summary>
        private readonly double jIntegralRadiusOverElementSize;

        /// <summary>
        /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
        /// rather than the using the J-integral method to predict it.
        /// </summary>
        private readonly PropagationLogger knownPropagation;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;

        /// <summary>
        /// The approximate ratio of the size of an element in the coarse region of the mesh to the size of an element in the 
        /// fine region of the mesh. Only applicable if <see cref="uniformMesh"/> == true.
        /// </summary>
        private readonly double ratioCoarseToFineElementSize;

        /// <summary>
        /// True to use same sized elements. False to refine the mesh near the crack path that is known roughly from previous 
        /// analysis.
        /// </summary>
        private readonly bool uniformMesh;

        /// <summary>
        /// If true, LSM will be used to describe the crack. If false, an explicit crack description (polyline) will be used.
        /// </summary>
        private readonly bool useLSM;

        private IMesh2D<XNode2D, XContinuumElement2D> mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="elementSize">The approximate size of the elements. If the mesh is not uniform, this applies to the 
        ///     finest region.</param>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private DCB(double elementSize, double growthLength, double jIntegralRadiusOverElementSize, 
            PropagationLogger knownPropagation, int maxIterations, double ratioCoarseToFineElementSize,
            bool uniformMesh, bool useLSM)
        {
            this.elementSize = elementSize;
            this.growthLength = growthLength;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.knownPropagation = knownPropagation;
            this.maxIterations = maxIterations;
            this.ratioCoarseToFineElementSize = ratioCoarseToFineElementSize;
            this.uniformMesh = uniformMesh;
            this.useLSM = useLSM;
        }
        
        /// <summary>
        /// The crack geometry description
        /// </summary>
        public IExteriorCrack Crack { get; private set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public Model2D Model { get; private set; }


        public IReadOnlyList<ICartesianPoint2D> Analyze(ISolver solver)
        {
            var actualPropagator = new Propagator(mesh, Crack, CrackTipPosition.Single, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            IPropagator propagator;
            if (knownPropagation != null) propagator = new FixedPropagator(knownPropagation, actualPropagator);
            else propagator = actualPropagator;
            var analysis = new QuasiStaticAnalysis(Model, mesh, Crack, solver, propagator, fractureToughness, maxIterations);
            return analysis.Analyze();
        }

        public void InitializeModel()
        {
            Model = new Model2D();
            CreateMesh();
            ApplyBoundaryConditions();
            InitializeCrack();
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(Model, 1e-6);

            // Fixed dofs
            foreach (var node in finder.FindNodesWithX(L))
            {
                Model.AddConstraint(node, DisplacementDOF.X, 0.0);
                Model.AddConstraint(node, DisplacementDOF.Y, 0.0);
            }

            // Loads
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            XNode2D topLeftNode = finder.FindNodeWith(0.0, h);
            Model.AddNodalLoad(bottomLeftNode, DisplacementDOF.Y, -load);
            Model.AddNodalLoad(topLeftNode, DisplacementDOF.Y, load);
        }

        private void CreateMesh()
        {
            // Mesh generation
            XNode2D[] nodes;
            List<XNode2D[]> elementConnectivity;
            if (uniformMesh) (nodes, elementConnectivity) = CreateUniformMesh();
            else (nodes, elementConnectivity) = CreateRectilinearMesh();

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
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
                Model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }

            // Mesh usable for crack-mesh interaction
            var boundary = new RectangularBoundary(0.0, L, 0.0, h);
            mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(Model.Nodes, Model.Elements, boundary);
        }

        private (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateUniformMesh()
        {
            int elementsPerXAxis = (int)(L / elementSize) + 1;
            int elementsPerYAxis = (int)(h / elementSize) + 1;
            var meshGenerator = new UniformRectilinearMeshGenerator(L, h, elementsPerXAxis, elementsPerYAxis);
            return meshGenerator.CreateMesh();
        }

        private (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateRectilinearMesh()
        {
            double fineElementSize = elementSize;
            double coarseElementSize = elementSize * ratioCoarseToFineElementSize;
            double[,] meshSizeAlongX = new double[,] 
            {
                { 0.0, coarseElementSize },
                { a, fineElementSize},
                //{ 3.0 / 5.0* L, fineElementSize},
                { 6.4, fineElementSize},
                { L, coarseElementSize}
            };
            double[,] meshSizeAlongY = new double[,] 
            {
                { 0.0, fineElementSize },
                //{ 0.0, fineElementSize},
                //{ 0.5, fineElementSize},
                { h/2 + h/10, elementSize * ratioCoarseToFineElementSize},
                { h, coarseElementSize}
            };

            var meshGenerator = new RectilinearMeshGenerator(meshSizeAlongX, meshSizeAlongY);
            return meshGenerator.CreateMesh();
        }

        private void InitializeCrack()
        {
            var crackVertex0 = new CartesianPoint2D(0.0, h / 2);
            var crackVertex1 = new CartesianPoint2D(a, h / 2);

            if (useLSM)
            {
                var lsmCrack = new BasicCrackLSM();
                lsmCrack.Mesh = mesh;

                // Create enrichments          
                lsmCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(lsmCrack);
                lsmCrack.CrackTipEnrichments = new CrackTipEnrichments2D(lsmCrack, CrackTipPosition.Single);

                // Mesh geometry interaction
                lsmCrack.InitializeGeometry(crackVertex0, crackVertex1);
                lsmCrack.UpdateGeometry(-dTheta, da);

                this.Crack = lsmCrack;
            }
            else
            {
                var explicitCrack = new BasicExplicitCrack2D();
                explicitCrack.Mesh = mesh;

                // Create enrichments          
                explicitCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(explicitCrack);
                explicitCrack.CrackTipEnrichments = new CrackTipEnrichments2D(explicitCrack, CrackTipPosition.Single);

                // Mesh geometry interaction
                explicitCrack.InitializeGeometry(crackVertex0, crackVertex1);
                explicitCrack.UpdateGeometry(-dTheta, da);

                this.Crack = explicitCrack;
            }
        }

        public class Builder
        {
            private readonly double elementSize;
            private readonly double growthLength;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="elementSize">The approximate size of the elements. If the mesh is not uniform, this applies to the 
            ///     finest region.</param>
            /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
            public Builder(double elementSize, double growthLength)
            {
                this.elementSize = elementSize;
                this.growthLength = growthLength;
            }

            /// <summary>
            /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
            /// radius is larger than the length of the crack segments.
            /// </summary>
            public double JintegralRadiusOverElementSize { get; set; } = 2.0;

            /// <summary>
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            public PropagationLogger KnownPropagation { get; set; } = null;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = int.MaxValue;

            /// <summary>
            /// The approximate ratio of the size of an element in the coarse region of the mesh to the size of an element in the 
            /// fine region of the mesh. Only applicable if <see cref="UniformMesh"/> == true.
            /// </summary>
            public double RatioCoarseToFineElementSize { get; set; } = 10;

            /// <summary>
            /// True to use same sized elements. False to refine the mesh near the crack path that is known roughly from previous 
            /// analysis.
            /// </summary>
            public bool UniformMesh { get; set; } = false;

            /// <summary>
            /// If true, LSM will be used to describe the crack. If false, an explicit crack description (polyline) will be used.
            /// </summary>
            public bool UseLSM { get; set; } = true;

            public DCB BuildBenchmark()
            {
                return new DCB(elementSize, growthLength, JintegralRadiusOverElementSize, KnownPropagation, MaxIterations,
                    RatioCoarseToFineElementSize, UniformMesh, UseLSM);
            }
        }
    }
}
