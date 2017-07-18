using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Analysis;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class DCBvariable
    {
        private static readonly double DIM_X = 60, DIM_Y = 20;
        private static readonly double E = 2e6, v = 0.3;
        private static readonly ICartesianPoint2D CRACK_MOUTH = new CartesianPoint2D(DIM_X, 0.5 * DIM_Y);

        public static void Main()
        {
            var test = new DCBvariable(15, 1);
            //test.CheckJintegralCountour();
            IVector solution = test.Solve();
            test.CheckSolution(solution);
            test.Propagate(solution);
        }


        private readonly SubmatrixChecker checker;

        //private readonly double elementSize;
        private readonly double jIntegralRadiusOverElementSize;

        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private Model2D model;
        private BasicExplicitCrack2D crack;

        private XNode2D bottomLeftNode;
        private XNode2D topLeftNode;
        private XNode2D bottomRightNode;
        private XNode2D topRightNode;

        private Propagator propagator;

        public DCBvariable(int elementsPerY, int jIntegralRadiusOverElementSize)
        {
            checker = new SubmatrixChecker(1e-4);

            //elementSize = DIM_Y / elementsPerY;
            //this.jIntegralRadius= jIntegralRadiusOverElementSize * this.elementSize;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;

            CreateModel(3 * elementsPerY, elementsPerY);
            HandleCrack();
        }

        private void CreateModel(int elementsPerX, int elementsPerY)
        {
            model = new Model2D();
            
            // Mesh
            var meshGenerator = new RectangularMeshGenerator(DIM_X, DIM_Y, elementsPerX, elementsPerY);
            Tuple<XNode2D[], List<XNode2D[]>> meshEntities = meshGenerator.CreateMesh();

            // Nodes
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);

            // Elements
            var integration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order4x4,
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4));

            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration, jIntegration));
            }

            // Boundary conditions
            var finder = new EntityFinder(model, 1e-6);
            bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            topLeftNode = finder.FindNodeWith(0.0, DIM_Y);
            bottomRightNode = finder.FindNodeWith(DIM_X, 0.0);
            topRightNode = finder.FindNodeWith(DIM_X, DIM_Y);

            model.AddConstraint(bottomLeftNode, StandardDOFType.X, 0.0);
            model.AddConstraint(bottomLeftNode, StandardDOFType.Y, 0.0);
            model.AddConstraint(topLeftNode, StandardDOFType.X, 0.0);
            model.AddConstraint(topLeftNode, StandardDOFType.Y, 0.0);
            model.AddConstraint(bottomRightNode, StandardDOFType.Y, -0.05);
            model.AddConstraint(topRightNode, StandardDOFType.Y, 0.05);
        }

        private void HandleCrack()
        {
            crack = new BasicExplicitCrack2D();
            crack.Mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements);

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack, new SignFunctionOpposite2D());
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack);

            // Mesh geometry interaction
            var crackTip = new CartesianPoint2D(0.5 * DIM_X, 0.5 * DIM_Y);
            crack.InitializeGeometry(CRACK_MOUTH, crackTip);
            crack.UpdateEnrichments();

            
        }

        private IVector Solve()
        {
            model.EnumerateDofs();
            var analysis = new LinearStaticAnalysis(model);
            analysis.Solve();
            return analysis.Solution;
        }

        private void CheckSolution(IVector solution)
        {
            var finder = new EntityFinder(model, 1e-6);
            List<XContinuumElement2D> mouthElements = finder.FindElementsThatContains(CRACK_MOUTH);
            XNode2D mouthBottomNode = mouthElements[0].Nodes[1];
            XNode2D mouthTopNode = mouthElements[0].Nodes[2];

            double uBotX = solution[model.DofEnumerator.GetFreeDofOf(mouthBottomNode, StandardDOFType.X)];
            double uBotY = solution[model.DofEnumerator.GetFreeDofOf(mouthBottomNode, StandardDOFType.Y)];
            double uTopX = solution[model.DofEnumerator.GetFreeDofOf(mouthTopNode, StandardDOFType.X)];
            double uTopY = solution[model.DofEnumerator.GetFreeDofOf(mouthTopNode, StandardDOFType.Y)];
            double aBotX = solution[model.DofEnumerator.GetArtificialDofOf(mouthBottomNode, crack.CrackBodyEnrichment.DOFs[0])];
            double aBotY = solution[model.DofEnumerator.GetArtificialDofOf(mouthBottomNode, crack.CrackBodyEnrichment.DOFs[1])];
            double aTopX = solution[model.DofEnumerator.GetArtificialDofOf(mouthTopNode, crack.CrackBodyEnrichment.DOFs[0])];
            double aTopY = solution[model.DofEnumerator.GetArtificialDofOf(mouthTopNode, crack.CrackBodyEnrichment.DOFs[1])];

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

        private void Propagate(IVector solution)
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
            propagator = new Propagator(crack.Mesh, crack, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(5));

            double[] totalConstrainedDisplacements = model.CalculateConstrainedDisplacements();
            double[] totalFreeDisplacements = new double[solution.Length];
            solution.CopyTo(totalFreeDisplacements, 0);

            double growthAngle, growthIncrement;
            propagator.Propagate(model, totalFreeDisplacements, totalConstrainedDisplacements, 
                out growthAngle, out growthIncrement);
            double jIntegral = (Math.Pow(propagator.Logger.SIFsMode1[0], 2) +
                Math.Pow(propagator.Logger.SIFsMode2[0], 2))
                / globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;

            Console.WriteLine("Propagation results:");
            propagator.Logger.PrintAnalysisStep(0);
            Console.WriteLine("J-integral = " + jIntegral);
            Console.WriteLine();
        }

        private void CheckJintegralCountour()
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
            propagator = new Propagator(crack.Mesh, crack, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.05));

            double radius = propagator.ComputeRadiusOfJintegralOuterContour();
            Circle2D outerContour = new Circle2D(crack.CrackTip, radius);
            IReadOnlyList<XContinuumElement2D> intersectedElements =
                crack.Mesh.FindElementsIntersectedByCircle(outerContour, crack.TipElements[0]);
        }
    }
}
