using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.Custom;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
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
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    internal class DoubleCantileverBeam
    {
        public const double beamHeight = 20;
        public const double beamLength = 60;
        public const int subdomainID = 0;
        public const double E = 2E6, v = 0.3;

        private readonly Rectangular2DBoundary beamBoundary;
        private readonly IIntegrationStrategy2D<XContinuumElement2D> integration, jIntegration;
        private readonly HomogeneousElasticMaterial2D material;

        public DoubleCantileverBeam()
        {
            this.material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
            this.integration = new IntegrationForCrackPropagation2D(
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)),
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)));
            this.jIntegration =
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(4, 4));
            this.beamBoundary = new Rectangular2DBoundary(0, beamLength, 0, beamHeight);
        }

        public TrackingExteriorCrackLSM Crack { get; private set; }
        public XModel Model { get; private set; }
        public Propagator Propagator { get; private set; }

        public void Create3x1Model()
        {
            double elementSize = 20.0;
            var model = new XModel();
            model.Subdomains[subdomainID] = new XSubdomain(subdomainID);
            this.Model = model;

            //Nodes
            model.Nodes.Add(new XNode(0, 0.0, 0.0));
            model.Nodes.Add(new XNode(1, elementSize, 0.0));
            model.Nodes.Add(new XNode(2, elementSize, elementSize));
            model.Nodes.Add(new XNode(3, 0.0, elementSize));
            model.Nodes.Add(new XNode(4, 2 * elementSize, 0.0));
            model.Nodes.Add(new XNode(5, 3 * elementSize, 0.0));
            model.Nodes.Add(new XNode(6, 3 * elementSize, elementSize));
            model.Nodes.Add(new XNode(7, 2 * elementSize, elementSize));

            // Elements
            XNode[][] connectivity = new XNode[3][];
            connectivity[0] = new XNode[] { model.Nodes[0], model.Nodes[1], model.Nodes[2], model.Nodes[3] };
            connectivity[1] = new XNode[] { model.Nodes[1], model.Nodes[4], model.Nodes[7], model.Nodes[2] };
            connectivity[2] = new XNode[] { model.Nodes[4], model.Nodes[5], model.Nodes[6], model.Nodes[7] };

            var factory = new XContinuumElement2DFactory(integration, jIntegration, material);
            var cells = new XContinuumElement2D[3];
            for (int e = 0; e < 3; ++e)
            {
                XContinuumElement2D element = factory.CreateElement(e, CellType.Quad4, connectivity[e]);
                cells[e] = element;
                model.Elements.Add(element);
                model.Subdomains[subdomainID].Elements.Add(model.Elements[e]);
            }

            // Mesh
            var mesh = new BidirectionalMesh2D<XNode, XContinuumElement2D>(model.Nodes,
                model.Elements.Select(e => (XContinuumElement2D)e).ToArray(), beamBoundary);

            // Boundary conditions
            model.Nodes[0].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            model.Nodes[0].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
            model.Nodes[3].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            model.Nodes[3].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
            model.Nodes[5].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = -0.05 });
            model.Nodes[6].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = +0.05 });

            // Create crack
            CreateCrack(elementSize, model, mesh, 2.0);
        }

        public void CreateModel(int numElementsX, int numElementsY, double jIntegralRadiusRatio)
        {
            var model = new XModel();
            model.Subdomains[subdomainID] = new XSubdomain(subdomainID);
            this.Model = model;

            // Mesh generation
            var meshGenerator = new UniformMeshGenerator2D<XNode>(0.0, 0.0, beamLength, beamHeight, numElementsX, numElementsY);
            (IReadOnlyList<XNode> nodes, IReadOnlyList<CellConnectivity<XNode>> elementConnectivities) = 
                meshGenerator.CreateMesh((id, x, y, z) => new XNode(id, x, y, z));

            // Nodes
            foreach (XNode node in nodes) model.Nodes.Add(node);

            // Elements
            var factory = new XContinuumElement2DFactory(integration, jIntegration, material);
            var cells = new XContinuumElement2D[elementConnectivities.Count];
            for (int e = 0; e < elementConnectivities.Count; ++e)
            {
                XContinuumElement2D element = factory.CreateElement(e, CellType.Quad4, elementConnectivities[e].Vertices);
                cells[e] = element;
                model.Elements.Add(element);
                model.Subdomains[subdomainID].Elements.Add(model.Elements[e]);
            }

            // Mesh
            var mesh = new BidirectionalMesh2D<XNode, XContinuumElement2D>(model.Nodes,
                model.Elements.Select(e => (XContinuumElement2D)e).ToArray(), beamBoundary);

            // Boundary conditions
            double tol = 1E-6;
            double L = DoubleCantileverBeam.beamLength;
            double H = DoubleCantileverBeam.beamHeight;
            XNode topRight = model.Nodes.Where(n => Math.Abs(n.X - L) <= tol && Math.Abs(n.Y - H) <= tol).First();
            XNode bottomRight = model.Nodes.Where(n => Math.Abs(n.X - L) <= tol && Math.Abs(n.Y) <= tol).First();
            topRight.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = +0.05 });
            bottomRight.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = -0.05 });
            foreach (XNode node in model.Nodes.Where(n => Math.Abs(n.X) <= tol))
            {
                node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
                node.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
            }

            // Create crack
            double elementSizeX = beamLength / numElementsX;
            double elementSizeY = beamHeight / numElementsY;
            CreateCrack(Math.Max(elementSizeX, elementSizeY), model, mesh, jIntegralRadiusRatio);
        }

        public (double jIintegral, double sifMode1) Propagate(Dictionary<int, Vector> freeDisplacements)
        {
            Crack.Propagate(freeDisplacements);
            double jIntegral = (Math.Pow(Propagator.Logger.SIFsMode1[0], 2) +
                Math.Pow(Propagator.Logger.SIFsMode2[0], 2))
                / material.HomogeneousEquivalentYoungModulus;

            return (jIntegral, Propagator.Logger.SIFsMode1[0]);
        }

        public (IVectorView globalU, IMatrixView globalK) SolveModel()
        {
            // Solver
            SkylineSolver solver = new SkylineSolver.Builder().BuildSolver(Model);
            solver.PreventFromOverwrittingSystemMatrices(); // Necessary to extract the stiffness matrix.

            // Problem type
            var provider = new ProblemStructural(Model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer(Model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(Model, solver, provider, childAnalyzer);

            // Run the anlaysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return (solver.LinearSystems[0].Solution, solver.LinearSystems[0].Matrix);
        }

        private void CreateCrack(double elementSize, XModel model, BidirectionalMesh2D<XNode, XContinuumElement2D> mesh, 
            double jIntegralRadiusRatio)
        {
            var propagator = new Propagator(mesh, jIntegralRadiusRatio, new HomogeneousMaterialAuxiliaryStates(material),
                new HomogeneousSIFCalculator(material), new MaximumCircumferentialTensileStressCriterion(),
                new ConstantIncrement2D(1.5 * jIntegralRadiusRatio * elementSize));
            var crack = new TrackingExteriorCrackLSM(propagator);
            crack.Mesh = mesh;

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack, new SignFunctionOpposite2D());
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);
            //crackTip = new CrackTip2D(CrackTip2D.TipCurvePosition.CurveStart, polyline, new SingleElementEnrichment(),
            //    2.0, new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //    new HomogeneousSIFCalculator(globalHomogeneousMaterial));

            // Mesh geometry interaction
            var crackMouth = new CartesianPoint(beamLength, 0.5 * beamHeight);
            var crackTip = new CartesianPoint(0.5 * beamLength, 0.5 * beamHeight);
            crack.InitializeGeometry(crackMouth, crackTip);
            crack.UpdateEnrichments(); //TODO: elements are not enriched properly.

            this.Crack = crack;
            this.Propagator = propagator;
        }
    }
}
