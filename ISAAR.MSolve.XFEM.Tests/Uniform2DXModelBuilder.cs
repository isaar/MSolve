using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.Custom;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition
{
    public class Uniform2DXModelBuilder
    {
        public enum BoundaryRegion
        {
            LeftSide, RightSide, UpperSide, LowerSide, UpperLeftCorner, UpperRightCorner, LowerLeftCorner, LowerRightCorner
        }

        private const double minX = 0.0, minY = 0.0;
        private List<(BoundaryRegion region, StructuralDof dof, double displacement)> prescribedDisplacements;
        private List<(BoundaryRegion region, StructuralDof dof, double load)> prescribedLoads;

        public Uniform2DXModelBuilder()
        {
            prescribedDisplacements = new List<(BoundaryRegion region, StructuralDof dof, double displacement)>();
            prescribedLoads = new List<(BoundaryRegion region, StructuralDof dof, double load)>();
        }

        public double DomainLengthX { get; set; } = 1.0;
        public double DomainLengthY { get; set; } = 1.0;
        public int NumSubdomainsX { get; set; } = 1;
        public int NumSubdomainsY { get; set; } = 1;
        public int NumTotalElementsX { get; set; } = 1;
        public int NumTotalElementsY { get; set; } = 1;

        public IIntegrationStrategy2D<XContinuumElement2D> QuadratureForStiffness { get; } =
            new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)));

        public IIntegrationStrategy2D<XContinuumElement2D> QuadratureForJintegral { get; } =
            new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(4, 4));

        public double YoungModulus { get; set; } = 1.0;

        /// <summary>
        /// Layout: left to right, then bottom to top. Example for 3x2 subdomains:
        /// ----------------  
        /// | E3 | E4 | E5 |
        /// ----------------
        /// | E0 | E1 | E2 |
        /// ----------------
        /// <see cref="YoungModuliOfSubdomains"/> = {{ E0, E1, E2 },{ E3, E4, E5 }}
        /// </summary>
        public double[,] YoungModuliOfSubdomains { get; set; } = null;

        public (XModel model, BidirectionalMesh2D<XNode, XContinuumElement2D> mesh) BuildModel()
        {
            // Generate global mesh
            double dx = DomainLengthX / NumTotalElementsX;
            double dy = DomainLengthY / NumTotalElementsY;
            double meshTolerance = 1E-10 * Math.Min(dx, dy);
            var meshGenerator = new UniformMeshGenerator2D<XNode>(0, 0, DomainLengthX, DomainLengthY,
                NumTotalElementsX, NumTotalElementsY);
            (IReadOnlyList<XNode> vertices, IReadOnlyList<CellConnectivity<XNode>> cells) = 
                meshGenerator.CreateMesh((id, x, y, z) => new XNode(id: id, x: x, y:  y, z: z ));

            // Define subdomain boundaries
            int numTotalSubdomains = NumSubdomainsX * NumSubdomainsY;
            var boundaries = new Rectangle[numTotalSubdomains];
            double subdomainLengthX = DomainLengthX / NumSubdomainsX;
            double subdomainLengthY = DomainLengthY / NumSubdomainsY;
            for (int j = 0; j < NumSubdomainsY; ++j)
            {
                double minY = j * subdomainLengthY;
                double maxY = (j + 1) * subdomainLengthY;
                for (int i = 0; i < NumSubdomainsX; ++i)
                {
                    double minX = i * subdomainLengthX;
                    double maxX = (i + 1) * subdomainLengthX;
                    boundaries[j * NumSubdomainsX + i] = new Rectangle(minX, minY, maxX, maxY);
                }
            }

            // Materials
            var youngModuli = new double[numTotalSubdomains];
            if (YoungModuliOfSubdomains == null)
            {
                for (int s = 0; s < numTotalSubdomains; ++s) youngModuli[s] = YoungModulus;
            }
            else
            {
                Debug.Assert(YoungModuliOfSubdomains.GetLength(0) == NumSubdomainsY
                    && YoungModuliOfSubdomains.GetLength(1) == NumSubdomainsX, "Materials do not match the subdomain layout");
                for (int j = 0; j < NumSubdomainsY; ++j)
                {
                    for (int i = 0; i < NumSubdomainsX; ++i)
                    {
                        youngModuli[j * NumSubdomainsX + i] = YoungModuliOfSubdomains[j, i];
                    }
                }
            }
            double thickness = 1.0;
            IMaterialField2D[] materials = youngModuli.Select(
                E => HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, 0.3, thickness)).ToArray();

            // Define model, subdomains, nodes
            var model = new XModel();
            for (int s = 0; s < numTotalSubdomains; ++s) model.Subdomains.Add(s, new XSubdomain(s));
            for (int n = 0; n < vertices.Count; ++n) model.Nodes.Add(vertices[n]);

            // Elements
            XContinuumElement2DFactory[] elementFactories = materials.Select(
                material => new XContinuumElement2DFactory(QuadratureForStiffness, QuadratureForJintegral, material)).ToArray();
            var elements = new XContinuumElement2D[cells.Count];
            for (int e = 0; e < cells.Count; ++e)
            {
                CellConnectivity<XNode> cell = cells[e];
                int NumSubdomainsContainingThis = 0;
                for (int s = 0; s < numTotalSubdomains; ++s)
                {
                    if (boundaries[s].Contains(cell, meshTolerance))
                    {
                        ++NumSubdomainsContainingThis;

                        // Create the element
                        XContinuumElement2D element = elementFactories[s].CreateElement(e, cell.CellType, cell.Vertices);
                        elements[e] = element;
                        model.Elements.Add(element);
                        model.Subdomains[s].Elements.Add(element);
                    }
                }
                Debug.Assert(NumSubdomainsContainingThis == 1);
            }

            // Mesh usable for crack-mesh interaction
            model.Boundary = new Rectangular2DBoundary(0.0, DomainLengthX, 0.0, DomainLengthY);
            var mesh = new BidirectionalMesh2D<XNode, XContinuumElement2D>(vertices, elements, model.Boundary);

            // Apply prescribed displacements
            foreach ((BoundaryRegion region, IDofType dof, double displacement) in prescribedDisplacements)
            {
                XNode[] nodes = FindBoundaryNodes(region, model, meshTolerance);
                foreach (XNode node in nodes) node.Constraints.Add(new Constraint() { DOF = dof, Amount = displacement });
            }

            // Apply prescribed loads
            foreach ((BoundaryRegion region, StructuralDof dof, double totalLoad) in prescribedLoads)
            {
                XNode[] nodes = FindBoundaryNodes(region, model, meshTolerance);
                double load = totalLoad / nodes.Length;
                foreach (XNode node in nodes) model.NodalLoads.Add(new NodalLoad(node, dof, load));
            }


            return (model, mesh);
        }

        /// <summary>
        /// </summary>
        /// <param name="load">Will be distributed evenly.</param>
        public void DistributeLoadAtNodes(BoundaryRegion region, StructuralDof dof, double load)
            => prescribedLoads.Add((region, dof, load));

        public void PrescribeDisplacement(BoundaryRegion region, StructuralDof dof, double displacement)
            => prescribedDisplacements.Add((region, dof, displacement));

        private XNode[] FindBoundaryNodes(BoundaryRegion region, XModel model, double tol)
        {
            double minX = 0.0, minY = 0.0, maxX = DomainLengthX, maxY = DomainLengthY; // for brevity

            IEnumerable<XNode> nodes;
            if (region == BoundaryRegion.LeftSide) nodes = model.Nodes.Where(node => Math.Abs(node.X - minX) <= tol);
            else if (region == BoundaryRegion.RightSide) nodes = model.Nodes.Where(node => Math.Abs(node.X - maxX) <= tol);
            else if (region == BoundaryRegion.LowerSide) nodes = model.Nodes.Where(node => Math.Abs(node.Y - minY) <= tol);
            else if (region == BoundaryRegion.UpperSide) nodes = model.Nodes.Where(node => Math.Abs(node.Y - maxY) <= tol);
            else if (region == BoundaryRegion.LowerLeftCorner)
            {
                nodes = model.Nodes.Where(node => (Math.Abs(node.X - minX) <= tol) && (Math.Abs(node.Y - minY) <= tol));
            }
            else if (region == BoundaryRegion.LowerRightCorner)
            {
                nodes = model.Nodes.Where(node => (Math.Abs(node.X - maxX) <= tol) && (Math.Abs(node.Y - minY) <= tol));
            }
            else if (region == BoundaryRegion.UpperLeftCorner)
            {
                nodes = model.Nodes.Where(node => (Math.Abs(node.X - minX) <= tol) && (Math.Abs(node.Y - maxY) <= tol));
            }
            else if (region == BoundaryRegion.UpperRightCorner)
            {
                nodes = model.Nodes.Where(node => (Math.Abs(node.X - maxX) <= tol) && (Math.Abs(node.Y - maxY) <= tol));
            }
            else throw new Exception("Should not have reached this code");

            return nodes.ToArray();
        }

        private class Rectangle
        {
            private readonly double minX, minY, maxX, maxY;

            internal Rectangle(double minX, double minY, double maxX, double maxY)
            {
                this.minX = minX;
                this.minY = minY;
                this.maxX = maxX;
                this.maxY = maxY;
            }

            public bool Contains(CellConnectivity<XNode> cell, double tol)
            {
                return cell.Vertices.All(node =>
                     (node.X >= minX - tol) && (node.X <= maxX + tol) && (node.Y >= minY - tol) && (node.Y <= maxY + tol));
            }
        }
    }
}
