using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class DCBpropagationBelytschko
    {
        private static readonly double L = 11.82, h = 3.94, t = 1.0; // in
        private static readonly double v = 0.3, E = 3e7; // psi=lbs/in^2
        private static readonly double load = 197; // lbs
        private static readonly double a = 3.95, da = 0.5; // in 
        private static readonly double dTheta = 5.71 * Math.PI / 180; // initial crack angle

        public static void Main()
        {
            // Parameters
            double jIntegralRadiusOverElementSize = 2.0;
            double fractureToughness = double.MaxValue;
            int maxIterations = 24;
            int defaultElementsPerY = 49;
            double defaultGrowthLength = 0.3;
            bool useLSM = true;
            //ParametricCrackLength(jIntegralRadiusOverElementSize, fractureToughness, maxIterations, defaultElementsPerY, useLSM);
            //ParametricMesh(jIntegralRadiusOverElementSize, fractureToughness, maxIterations, defaultGrowthLength, useLSM);
            GridSearch(jIntegralRadiusOverElementSize, fractureToughness, maxIterations, useLSM);
        }

        private static void ParametricCrackLength(double jIntegralRadiusOverElementSize, double fractureToughness, 
            int maxIterations, int elementsPerY, bool useLSM)
        {
            double[] growthLengths = new double[] { 0.1, 0.2, 0.3, 0.4, 0.5 };
            Console.WriteLine("------------------------------------ Parametric Growth Length ------------------------------------");
            for (int i = 0; i < growthLengths.Length; ++i)
            {
                Console.WriteLine("------------------ Mesh: {0}x{1} , Growth length = {2} ------------------",
                    elementsPerY, 3 * elementsPerY, growthLengths[i]);
                var benchmark = new DCBpropagationBelytschko(elementsPerY, jIntegralRadiusOverElementSize,
                    fractureToughness, maxIterations, growthLengths[i], useLSM);
                try
                {
                    IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze();
                    Console.WriteLine("Crack path:");
                    foreach (var point in crackPath)
                    {
                        Console.WriteLine("{0} {1}", point.X, point.Y);
                    }
                }
                catch (Exception e)
                {
                    //Console.WriteLine(e.Message);
                    throw e;
                }
                Console.WriteLine();
            }
        }

        private static void ParametricMesh(double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations, double propagationLength, bool useLSM)
        {
            int[] elementsPerY = new int[] { 15, 23, 37, 51 };
            Console.WriteLine("------------------------------------ Parametric Mesh ------------------------------------");
            for (int i = 0; i < elementsPerY.Length; ++i)
            {
                Console.WriteLine("------------------ Mesh: {0}x{1} , Growth length = {2} ------------------",
                    elementsPerY[i], 3 * elementsPerY[i], propagationLength);
                var benchmark = new DCBpropagationBelytschko(elementsPerY[i], jIntegralRadiusOverElementSize,
                    fractureToughness, maxIterations, propagationLength, useLSM);
                try
                {
                    IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze();
                    Console.WriteLine("Crack path:");
                    foreach (var point in crackPath)
                    {
                        Console.WriteLine("{0} {1}", point.X, point.Y);
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
                Console.WriteLine();
            }
        }

        private static void GridSearch(double jIntegralRadiusOverElementSize, double fractureToughness, 
            int maxIterations, bool useLSM)
        {
            for (int elems = 15; elems <= 60; ++elems)
            {
                //double[] growthLengths = new double[] { 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 };
                double[] growthLengths = new double[] { 0.2 };
                for (int i = 0; i < growthLengths.Length; ++i)
                {
                    Console.WriteLine("------------------ Mesh: {0}x{1}={2} , Growth length = {3} ------------------",
                        elems, 3 * elems, 3 * elems * elems, growthLengths[i]);
                    var benchmark = new DCBpropagationBelytschko(elems, jIntegralRadiusOverElementSize,
                        fractureToughness, maxIterations, growthLengths[i], useLSM);
                    try
                    {
                        IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze();
                        Console.WriteLine("Crack path:");
                        foreach (var point in crackPath)
                        {
                            Console.WriteLine("{0} {1}", point.X, point.Y);
                        }
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine(e.Message);
                    }
                    Console.WriteLine();
                }
            }
        }

        private Model2D model;
        private readonly bool useLSM;
        private IExteriorCrack crack;
        private IMesh2D<XNode2D, XContinuumElement2D> mesh;
        private IIntegrationStrategy2D<XContinuumElement2D> integration;
        private IIntegrationStrategy2D<XContinuumElement2D> jIntegration;
        private readonly double jIntegralRadiusOverElementSize;
        private Propagator propagator;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private readonly double fractureToughness;
        private readonly int maxIterations;
        private readonly double growthLength;

        public DCBpropagationBelytschko(int elementsPerY, double jIntegralRadiusOverElementSize,
            double fractureToughness, int maxIterations, double growthLength, bool useLSM)
        {
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
            this.growthLength = growthLength;
            this.useLSM = useLSM;

            CreateModel(elementsPerY);
        }

        private void CreateModel(int elementsPerY)
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
            crack = new BasicExplicitCrack2D();
            model = new Model2D();
            HandleIntegrations();
            CreateMesh(elementsPerY);
            ApplyBoundaryConditions();
            HandleCrack();
        }

        private void HandleIntegrations()
        {
            integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);
        }

        private void CreateMesh(int elementsPerY)
        {
            var meshGenerator = new UniformRectilinearMeshGenerator(L, h, 3 * elementsPerY, elementsPerY);
            Tuple<XNode2D[], List<XNode2D[]>> meshEntities = meshGenerator.CreateMesh();
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);
            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }

            var boundary = new RectangularBoundary(0.0, L, 0.0, h);
            mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements, boundary);
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(model, 1e-6);

            // Fixed dofs
            foreach (var node in finder.FindNodesWithX(L))
            {
                model.AddConstraint(node, StandardDOFType.X, 0.0);
                model.AddConstraint(node, StandardDOFType.Y, 0.0);
            }

            // Loads
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            XNode2D topLeftNode = finder.FindNodeWith(0.0, h);
            model.AddNodalLoad(bottomLeftNode, StandardDOFType.Y, -load);
            model.AddNodalLoad(topLeftNode, StandardDOFType.Y, load);
        }

        private void HandleCrack()
        {
            if (useLSM)
            {
                var lsmCrack = new BasicCrackLSM();
                this.crack = lsmCrack;
                lsmCrack.Mesh = mesh;

                // Create enrichments          
                lsmCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(lsmCrack);
                lsmCrack.CrackTipEnrichments = new CrackTipEnrichments2D(lsmCrack, CrackTipPosition.Single);

                // Mesh geometry interaction
                var crackVertex0 = new CartesianPoint2D(0.0, h / 2);
                var crackVertex1 = new CartesianPoint2D(a, h / 2);
                lsmCrack.InitializeGeometry(crackVertex0, crackVertex1);
                lsmCrack.UpdateGeometry(-dTheta, da);
            }
            else
            {
                var explicitCrack = new BasicExplicitCrack2D();
                this.crack = explicitCrack;
                explicitCrack.Mesh = mesh;

                // Create enrichments          
                explicitCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(explicitCrack);
                explicitCrack.CrackTipEnrichments = new CrackTipEnrichments2D(explicitCrack, CrackTipPosition.Single);

                // Mesh geometry interaction
                var crackVertex0 = new CartesianPoint2D(0.0, h / 2);
                var crackVertex1 = new CartesianPoint2D(a, h / 2);
                explicitCrack.InitializeGeometry(crackVertex0, crackVertex1);
                explicitCrack.UpdateGeometry(-dTheta, da);
            }
            
        }

        private IReadOnlyList<ICartesianPoint2D> Analyze()
        {
            propagator = new Propagator(mesh, crack, CrackTipPosition.Single, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            QuasiStaticAnalysis analysis = new QuasiStaticAnalysis(model, mesh, crack,
                propagator, fractureToughness, maxIterations);
            return analysis.Analyze();
        }
    }
}
