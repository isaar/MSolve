using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.SuiteSparse;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;
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
        private static readonly bool useLSM = false;
        private static readonly bool uniformMesh = false;
        private static readonly double coarseOverFine = 10;

        public static void Main()
        {
            // Parameters
            double jIntegralRadiusOverElementSize = 2.0;
            double fractureToughness = double.MaxValue;
            int maxIterations = 100;
            SingleTest(jIntegralRadiusOverElementSize, fractureToughness, maxIterations);
            //ParametricCrackLength(jIntegralRadiusOverElementSize, fractureToughness, maxIterations);
            //ParametricMesh(jIntegralRadiusOverElementSize, fractureToughness, maxIterations);
            //GridSearch(jIntegralRadiusOverElementSize, fractureToughness, maxIterations);
        }

        private static void SingleTest(double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations)
        {
            double growthLength = 0.05;
            double fineElementSize = 0.08;

            var benchmark = new DCBpropagationBelytschko(fineElementSize, jIntegralRadiusOverElementSize,
                        fractureToughness, maxIterations, growthLength);
            Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                fineElementSize, benchmark.model.Elements.Count, growthLength);

            //Print
            //VTKWriter writer = new VTKWriter(benchmark.model);
            //writer.InitializeFile("dcb_transfinite");
            //writer.CloseCurrentFile();

            IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze();
            Console.WriteLine("Crack path:");
            foreach (var point in crackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
        }

        private static void ParametricCrackLength(double jIntegralRadiusOverElementSize, double fractureToughness, 
            int maxIterations)
        {
            double fineElementSize = 0.038;
            double[] growthLengths = new double[] { 0.05, 0.1, 0.2, 0.4 };
            Console.WriteLine("------------------------------------ Parametric Growth Length ------------------------------------");
            for (int i = 0; i < growthLengths.Length; ++i)
            {
                var benchmark = new DCBpropagationBelytschko(fineElementSize, jIntegralRadiusOverElementSize,
                    fractureToughness, maxIterations, growthLengths[i]);
                Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                    fineElementSize, benchmark.model.Elements.Count, growthLengths[i]);
                
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
                    //throw e;
                }
                Console.WriteLine();
            }
        }

        private static void ParametricMesh(double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations)
        {
            double propagationLength = 0.3;
            //double[] fineElementSizes = new double[] { 0.046 };
            double[] fineElementSizes = new double[] { 0.046, 0.1, 0.15, 0.2, 0.25 };
            Console.WriteLine("------------------------------------ Parametric Mesh ------------------------------------");
            for (int i = 0; i < fineElementSizes.Length; ++i)
            {
                var benchmark = new DCBpropagationBelytschko(fineElementSizes[i], jIntegralRadiusOverElementSize,
                    fractureToughness, maxIterations, propagationLength);
                Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                    fineElementSizes[i], benchmark.model.Elements.Count, propagationLength);

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
            int maxIterations)
        {
            double[] growthLengths = new double[] {0.2, 0.25, 0.3};
            double[] fineElementSizes = new double[] { 0.035, 0.036, 0.037, 0.038, 0.039, 0.040, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049 };
            //double[] growthLengths = new double[] { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
            //double[] fineElementSizes = new double[] { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 };
            //double[] fineElementSizes = new double[] {0.035, 0.036, 0.037, 0.038, 0.039, 0.040, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049 };

            for (int j = 0; j < growthLengths.Length; ++j)
            {
                for (int i = 0; i < fineElementSizes.Length; ++i)
                {
                    var benchmark = new DCBpropagationBelytschko(fineElementSizes[i], jIntegralRadiusOverElementSize,
                        fractureToughness, maxIterations, growthLengths[j]);
                    Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                        fineElementSizes[i], benchmark.model.Elements.Count, growthLengths[j]);
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

        public Model2D model;
        private TrackingExteriorCrackLSM crack;
        private BiMesh2D mesh;
        private IIntegrationStrategy2D<XContinuumElement2D> integration;
        private IIntegrationStrategy2D<XContinuumElement2D> jIntegration;
        private readonly double jIntegralRadiusOverElementSize;
        private Propagator propagator;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private readonly double fractureToughness;
        private readonly int maxIterations;
        private readonly double growthLength;

        public DCBpropagationBelytschko(double fineElementSize, double jIntegralRadiusOverElementSize,
            double fractureToughness, int maxIterations, double growthLength)
        {
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
            this.growthLength = growthLength;

            CreateModel(fineElementSize);
        }

        private void CreateModel(double fineElementSize)
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
            model = new Model2D();
            HandleIntegrations();
            CreateMesh(fineElementSize);
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

        private void CreateMesh(double fineElementSize)
        {
            XNode2D[] nodes;
            List<XNode2D[]> elementNodes;
            if (uniformMesh) (nodes, elementNodes) = CreateUniformMesh(fineElementSize);
            else (nodes, elementNodes) = CreateRectilinearMesh(fineElementSize);

            foreach (XNode2D node in nodes) model.AddNode(node);
            foreach (XNode2D[] element in elementNodes)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, element, materialField,
                    integration, jIntegration));
            }

            var boundary = new RectangularBoundary(0.0, L, 0.0, h);
            mesh = new BiMesh2D(model.Nodes, model.Elements, boundary);
        }


        private (XNode2D[] nodes, List<XNode2D[]> elementNodes) CreateUniformMesh(double fineElementSize)
        {
            int elementsPerXAxis = (int)(L / fineElementSize) + 1;
            int elementsPerYAxis = (int)(h / fineElementSize) + 1;
            var meshGenerator = new UniformRectilinearMeshGenerator(L, h, elementsPerXAxis, elementsPerYAxis);
            return meshGenerator.CreateMesh();
        }

        private (XNode2D[] nodes, List<XNode2D[]> elementNodes) CreateRectilinearMesh(double fineElementSize)
        {
            double[,] meshSizeAlongX = new double[,] {
                { 0.0, fineElementSize * coarseOverFine},
                { a, fineElementSize},
                //{ 3.0 / 5.0* L, fineElementSize},
                { 6.4, fineElementSize},
                { L, fineElementSize * coarseOverFine} };
            double[,] meshSizeAlongY = new double[,] {
                { 0.0, fineElementSize },
                //{ 0.0, fineElementSize * coarseOverFine},
                //{ 0.5, fineElementSize},
                { h/2 + h/10, fineElementSize},
                { h, fineElementSize * coarseOverFine} };

            var meshGenerator = new RectilinearMeshGenerator(meshSizeAlongX, meshSizeAlongY);
            return meshGenerator.CreateMesh();
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(model, 1e-6);

            // Fixed dofs
            foreach (var node in finder.FindNodesWithX(L))
            {
                model.AddConstraint(node, DisplacementDof.X, 0.0);
                model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }

            // Loads
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            XNode2D topLeftNode = finder.FindNodeWith(0.0, h);
            model.AddNodalLoad(bottomLeftNode, DisplacementDof.Y, -load);
            model.AddNodalLoad(topLeftNode, DisplacementDof.Y, load);
        }

        private void HandleCrack()
        {
            propagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
               new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
               new HomogeneousSIFCalculator(globalHomogeneousMaterial),
               new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            if (useLSM)
            {
                var lsmCrack = new TrackingExteriorCrackLSM(propagator);
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
                throw new NotImplementedException();
                var explicitCrack = new BasicExplicitCrack2D();
                //this.crack = explicitCrack;
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
            var solver = new SkylineSolver(model);
            QuasiStaticAnalysis analysis = new QuasiStaticAnalysis(model, mesh, crack, solver, fractureToughness, maxIterations);
            analysis.Analyze();
            return crack.CrackPath;
        }
    }
}
