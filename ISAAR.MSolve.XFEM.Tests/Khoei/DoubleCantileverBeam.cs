using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
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
        public const int subdomainID = 0;

        public DoubleCantileverBeam()
        {
        }

        public TrackingExteriorCrackLSM Crack { get; private set; }
        public XModel Model { get; private set; }

        public static DoubleCantileverBeam Create3x1() => new DoubleCantileverBeam();

        public void Create3x1Model()
        {
            double elementSize = 20.0;
            var model = new XModel();
            model.Subdomains[subdomainID] = new XSubdomain(subdomainID);

            //Nodes
            model.Nodes.Add(new XNode(0, 0.0, 0.0));
            model.Nodes.Add(new XNode(1, elementSize, 0.0));
            model.Nodes.Add(new XNode(2, elementSize, elementSize));
            model.Nodes.Add(new XNode(3, 0.0, elementSize));
            model.Nodes.Add(new XNode(4, 2 * elementSize, 0.0));
            model.Nodes.Add(new XNode(5, 3 * elementSize, 0.0));
            model.Nodes.Add(new XNode(6, 3 * elementSize, elementSize));
            model.Nodes.Add(new XNode(7, 2 * elementSize, elementSize));

            // Material
            double E = 2e6, v = 0.3;
            var commonMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

            // Integration rules
            IIntegrationStrategy2D<XContinuumElement2D> integration, jIntegration;
            integration = new IntegrationForCrackPropagation2D(
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)),
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)));
            jIntegration =
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(4, 4));

            // Elements
            XNode[][] connectivity = new XNode[3][];
            connectivity[0] = new XNode[] { model.Nodes[0], model.Nodes[1], model.Nodes[2], model.Nodes[3] };
            connectivity[1] = new XNode[] { model.Nodes[1], model.Nodes[4], model.Nodes[7], model.Nodes[2] };
            connectivity[2] = new XNode[] { model.Nodes[4], model.Nodes[5], model.Nodes[6], model.Nodes[7] };

            var factory = new XContinuumElement2DFactory(integration, jIntegration, commonMaterial);
            var cells = new XContinuumElement2D[3];
            for (int e = 0; e < 3; ++e)
            {
                XContinuumElement2D element = factory.CreateElement(e, CellType.Quad4, connectivity[e]);
                cells[e] = element;
                model.Elements.Add(element);
                model.Subdomains[subdomainID].Elements.Add(model.Elements[e]);
            }

            // Mesh
            var boundary = new Rectangular2DBoundary(0, 3 * elementSize, 0, elementSize);
            var mesh = new BidirectionalMesh2D<XNode, XContinuumElement2D>(model.Nodes,
                model.Elements.Select(e => (XContinuumElement2D)e).ToArray(), boundary);

            // Boundary conditions
            model.Nodes[0].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            model.Nodes[0].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
            model.Nodes[3].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX, Amount = 0.0 });
            model.Nodes[3].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = 0.0 });
            model.Nodes[5].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = -0.05 });
            model.Nodes[6].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY, Amount = +0.05 });

            this.Model = model;
            this.Crack = CreateCrack(elementSize, model, mesh, commonMaterial);
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

        private TrackingExteriorCrackLSM CreateCrack(double elementSize, XModel model,
            BidirectionalMesh2D<XNode, XContinuumElement2D> mesh, HomogeneousElasticMaterial2D material)
        {
            var propagator = new Propagator(mesh, 2.0, new HomogeneousMaterialAuxiliaryStates(material),
                new HomogeneousSIFCalculator(material), new MaximumCircumferentialTensileStressCriterion(),
                new ConstantIncrement2D(2.0 * elementSize));
            var crack = new TrackingExteriorCrackLSM(propagator);
            crack.Mesh = mesh;

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack, new SignFunctionOpposite2D());
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);
            //crackTip = new CrackTip2D(CrackTip2D.TipCurvePosition.CurveStart, polyline, new SingleElementEnrichment(),
            //    2.0, new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //    new HomogeneousSIFCalculator(globalHomogeneousMaterial));

            // Mesh geometry interaction
            var crackMouth = new CartesianPoint(3 * elementSize, 0.5 * elementSize);
            var crackTip = new CartesianPoint(1.5 * elementSize, 0.5 * elementSize);
            crack.InitializeGeometry(crackMouth, crackTip);
            crack.UpdateEnrichments(); //TODO: elements are not enriched properly.

            return crack;
        }
    }
}
