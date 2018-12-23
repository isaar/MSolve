using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.SuiteSparse;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class DCBvariable
    {
        private static readonly double DIM_X = 60, DIM_Y = 20;
        private static readonly double E = 2e6, v = 0.3;
        private static readonly ICartesianPoint2D CRACK_MOUTH = new CartesianPoint2D(DIM_X, 0.5 * DIM_Y);
        private static readonly bool structuredMesh = true;
        private static readonly string triMeshFile = "dcb_tri.msh";
        private static readonly string quadMeshFile = "dcb_quad.msh";
        private static readonly bool integrationWithTriangles = false;
        //private static readonly double triangleOverElementArea = double.PositiveInfinity; Can't refine yet.

        public static void Main()
        {
            // Error in 9x27 and higher for suitesparse solver, but not for all j-integral radii
            int[] meshElements = new int[] {  15, 25, 45 }; 
            double[] jIntegralRadiiOverElementSize = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var watch = new Stopwatch();
            watch.Start();
            Console.WriteLine("---------------------- Results ---------------------");
            for (int i = 0; i < meshElements.Length; ++i)
            {
                for (int j = 0; j < jIntegralRadiiOverElementSize.Length; ++j)
                {
                    var test = new DCBvariable(meshElements[i], jIntegralRadiiOverElementSize[j]);
                    //test.CheckJintegralCountour();
                    Vector solution = test.Solve();
                    //test.CheckSolution(solution);
                    Tuple<double, double> results = test.Propagate(solution);

                    Console.WriteLine("Mesh = ({0}x{1}), J-integral radius / element size = {2}: J = {3}, KI = {4}",
                    meshElements[i], 3 * meshElements[i], jIntegralRadiiOverElementSize[j], results.Item1, results.Item2);
                }
            }
            watch.Stop();
            double duration = watch.ElapsedMilliseconds / 1000.0;
            Console.WriteLine($"----- Time elapsed = {duration} s -----");
        }

        private readonly SubmatrixChecker checker;
        //private readonly double elementSize;
        private readonly double jIntegralRadiusOverElementSize;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private Model2D model;
        private IDofOrderer dofOrderer;
        private TrackingExteriorCrackLSM crack;
        private IIntegrationStrategy2D<XContinuumElement2D> integration;
        private IIntegrationStrategy2D<XContinuumElement2D> jIntegration;
        //private XNode2D bottomLeftNode;
        //private XNode2D topLeftNode;
        private XNode2D bottomRightNode;
        private XNode2D topRightNode;

        private Propagator propagator;

        public DCBvariable(int elementsPerY, double jIntegralRadiusOverElementSize)
        {
            checker = new SubmatrixChecker(1e-4);

            //elementSize = DIM_Y / elementsPerY;
            //this.jIntegralRadius= jIntegralRadiusOverElementSize * this.elementSize;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;

            CreateModel(elementsPerY);
            
        }

        private void CreateModel(int elementsPerY)
        {
            // TODO: Fix the construction order. The crack needed for integration schemes. However the mesh depends
            // on the elements than need the integration schemes. Perhaps pass the model to the crack and acces s
            // the mesh through the model.
            model = new Model2D();
            HandleIntegrations();
            if (structuredMesh) CreateStructuredMesh(elementsPerY);
            else CreateUnstructuredMesh();
            ApplyBoundaryConditions();
            HandleCrack();
        }

        private void HandleIntegrations()
        {
            if (integrationWithTriangles)
            {
                throw new NotImplementedException("The crack is not defined yet, thus cannot be used for creating triangles."); ; 
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

        private void CreateStructuredMesh(int elementsPerY)
        {
            var meshGenerator = new UniformRectilinearMeshGenerator(DIM_X, DIM_Y, 3 * elementsPerY, elementsPerY);
            (XNode2D[] nodes, List<XNode2D[]> elementNodes) meshEntities = meshGenerator.CreateMesh();
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);
            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }
        }

        private void CreateUnstructuredMesh()
        {
            //var meshReader = new GmshReader(triMeshFile);
            var meshReader = new GmshReader(quadMeshFile);
            Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> meshEntities = meshReader.ReadMesh();
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);
            foreach (var element in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                model.AddElement(new XContinuumElement2D(element.ElementType, element.Nodes, materialField, 
                    integration, jIntegration));
            }
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(model, 1e-6);

            // Fixed dofs
            foreach (var node in finder.FindNodesWithX(0.0))
            {
                model.AddConstraint(node, DisplacementDof.X, 0.0);
                model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }
            //bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            //topLeftNode = finder.FindNodeWith(0.0, DIM_Y);
            //model.AddConstraint(bottomLeftNode, DisplacementDof.X, 0.0);
            //model.AddConstraint(bottomLeftNode, DisplacementDof.Y, 0.0);
            //model.AddConstraint(topLeftNode, DisplacementDof.X, 0.0);
            //model.AddConstraint(topLeftNode, DisplacementDof.Y, 0.0);

            // Prescribed displacements
            bottomRightNode = finder.FindNodeWith(DIM_X, 0.0);
            topRightNode = finder.FindNodeWith(DIM_X, DIM_Y);
            model.AddConstraint(bottomRightNode, DisplacementDof.Y, -0.05);
            model.AddConstraint(topRightNode, DisplacementDof.Y, 0.05);
        }

        private void HandleCrack()
        {
            var boundary = new RectangularBoundary(0.0, DIM_X, 0.0, DIM_Y);
            var mesh = new BiMesh2D(model.Nodes, model.Elements, boundary);

            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
            propagator = new Propagator(crack.Mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(5));


            crack = new TrackingExteriorCrackLSM(propagator);
            crack.Mesh = mesh;
            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack, new SignFunctionOpposite2D());
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);

            // Mesh geometry interaction
            var crackTip = new CartesianPoint2D(0.5 * DIM_X, 0.5 * DIM_Y);
            crack.InitializeGeometry(CRACK_MOUTH, crackTip);
            crack.UpdateEnrichments();
        }

        private Vector Solve()
        {
            var solver = new CholeskySuiteSparseSolver(model);
            //var solver = new SkylineSolver(model);
            //var solver = new PCGSolver(model, 1, 1e-8);
            solver.Initialize();
            solver.Solve();
            dofOrderer = solver.DofOrderer;
            return (Vector)solver.Solution;
        }

        private void CheckSolution(Vector solution)
        {
            var finder = new EntityFinder(model, 1e-6);
            List<XContinuumElement2D> mouthElements = finder.FindElementsThatContains(CRACK_MOUTH);
            XNode2D mouthBottomNode = mouthElements[0].Nodes[1];
            XNode2D mouthTopNode = mouthElements[0].Nodes[2];

            double uBotX = solution[dofOrderer.GetStandardDofOf(mouthBottomNode, DisplacementDof.X)];
            double uBotY = solution[dofOrderer.GetStandardDofOf(mouthBottomNode, DisplacementDof.Y)];
            double uTopX = solution[dofOrderer.GetStandardDofOf(mouthTopNode, DisplacementDof.X)];
            double uTopY = solution[dofOrderer.GetStandardDofOf(mouthTopNode, DisplacementDof.Y)];
            double aBotX = solution[dofOrderer.GetEnrichedDofOf(mouthBottomNode, crack.CrackBodyEnrichment.Dofs[0])];
            double aBotY = solution[dofOrderer.GetEnrichedDofOf(mouthBottomNode, crack.CrackBodyEnrichment.Dofs[1])];
            double aTopX = solution[dofOrderer.GetEnrichedDofOf(mouthTopNode, crack.CrackBodyEnrichment.Dofs[0])];
            double aTopY = solution[dofOrderer.GetEnrichedDofOf(mouthTopNode, crack.CrackBodyEnrichment.Dofs[1])];

            Console.WriteLine("Solution results: For the element containing the crack mouth:");
            Console.WriteLine("Bottom right node, standard dof x: u = " + uBotX);
            Console.WriteLine("Bottom right node, standard dof y: u = " + uBotY);
            Console.WriteLine("Top right node, standard dof x: u = " + uTopX);
            Console.WriteLine("Top right node, standard dof y: u = " + uTopY);
            Console.WriteLine("Bottom right node, enriched dof x: u = " + aBotX);
            Console.WriteLine("Bottom right node, enriched dof y: u = " + aBotY);
            Console.WriteLine("Top right node, enriched dof x: u = " + aTopX);
            Console.WriteLine("Top right node, enriched dof y: u = " + aTopY);
            Console.WriteLine();
        }

        private Tuple<double, double> Propagate(Vector solution)
        {
            Vector totalConstrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);

            crack.Propagate(dofOrderer, solution, totalConstrainedDisplacements);
            double jIntegral = (Math.Pow(propagator.Logger.SIFsMode1[0], 2) +
                Math.Pow(propagator.Logger.SIFsMode2[0], 2))
                / globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;

            return new Tuple<double, double>(jIntegral, propagator.Logger.SIFsMode1[0]);

            //Console.WriteLine("Propagation results:");
            //propagator.Logger.PrintAnalysisStep(0);
            //Console.WriteLine("J-integral = " + jIntegral);
            //Console.WriteLine();
        }

        private void CheckJintegralCountour()
        {
            throw new NotImplementedException();
            //globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
            //propagator = new Propagator(crack.Mesh, crack, CrackTipPosition.Single, jIntegralRadiusOverElementSize,
            //    new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //    new HomogeneousSIFCalculator(globalHomogeneousMaterial),
            //    new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.05));

            //double radius = propagator.ComputeRadiusOfJintegralOuterContour();
            //Circle2D outerContour = new Circle2D(crack.GetCrackTip(CrackTipPosition.Single), radius);
            //IReadOnlyList<XContinuumElement2D> intersectedElements =
            //    crack.Mesh.FindElementsIntersectedByCircle(outerContour, 
            //    crack.GetTipElements(CrackTipPosition.Single)[0]);
        }
    }
}
