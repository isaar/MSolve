using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition
{
    public class Uniform2DModelBuilder
    {
        public enum BoundaryRegion
        {
            LeftSide, RightSide, UpperSide, LowerSide, UpperLeftCorner, UpperRightCorner, LowerLeftCorner, LowerRightCorner
        }

        private const double minX = 0.0, minY = 0.0;
        private List<(BoundaryRegion region, DOFType dof, double displacement)> prescribedDisplacements;
        private List<(BoundaryRegion region, DOFType dof, double load)> prescribedLoads;

        public Uniform2DModelBuilder()
        {
            prescribedDisplacements = new List<(BoundaryRegion region, DOFType dof, double displacement)>();
            prescribedLoads = new List<(BoundaryRegion region, DOFType dof, double load)>();
        }

        public double DomainLengthX { get; set; } = 1.0;
        public double DomainLengthY { get; set; } = 1.0;
        public int NumSubdomainsX { get; set; } = 1;
        public int NumSubdomainsY { get; set; } = 1;
        public int NumTotalElementsX { get; set; } = 1;
        public int NumTotalElementsY { get; set; } = 1;

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

        public Model_v2 BuildModel()
        {
            // Generate global mesh
            double dx = DomainLengthX / NumTotalElementsX;
            double dy = DomainLengthY / NumTotalElementsY;
            double meshTolerance = 1E-10 * Math.Min(dx, dy);
            var meshGenerator = new UniformMeshGenerator2D_v2(0, 0, DomainLengthX, DomainLengthY,
                NumTotalElementsX, NumTotalElementsY);
            (IReadOnlyList<Node_v2> vertices, IReadOnlyList<CellConnectivity_v2> cells) = meshGenerator.CreateMesh();

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
            var dynamicProperties = new DynamicMaterial(1.0, 0.0, 0.0);
            ElasticMaterial2D_v2[] materials = youngModuli.Select(
                E => new ElasticMaterial2D_v2(StressState2D.PlaneStress) { YoungModulus = E, PoissonRatio = 0.3 }).ToArray();

            // Define model, subdomains, nodes
            var model = new Model_v2();
            for (int s = 0; s < numTotalSubdomains; ++s) model.SubdomainsDictionary.Add(s, new Subdomain_v2(s));
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Elements
            ContinuumElement2DFactory[] elementFactories = materials.Select(
                material => new ContinuumElement2DFactory(thickness, material, dynamicProperties)).ToArray();
            for (int e = 0; e < cells.Count; ++e)
            {
                int NumSubdomainsContainingThis = 0;
                for (int s = 0; s < numTotalSubdomains; ++s)
                {
                    if (boundaries[s].Contains(cells[e], meshTolerance))
                    {
                        ++NumSubdomainsContainingThis;

                        // Create the element
                        ContinuumElement2D element = elementFactories[s].CreateElement(cells[e].CellType, cells[e].Vertices);
                        var elementWrapper = new Element_v2() { ID = e, ElementType = element };
                        foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                        model.ElementsDictionary.Add(e, elementWrapper);
                        model.SubdomainsDictionary[s].Elements.Add(elementWrapper);
                    }
                }
                Debug.Assert(NumSubdomainsContainingThis == 1);
            }

            // Apply prescribed displacements
            foreach ((BoundaryRegion region, DOFType dof, double displacement) in prescribedDisplacements)
            {
                Node_v2[] nodes = FindBoundaryNodes(region, model, meshTolerance);
                foreach (Node_v2 node in nodes) node.Constraints.Add(new Constraint() { DOF = dof, Amount = displacement });
            }

            // Apply prescribed loads
            foreach ((BoundaryRegion region, DOFType dof, double totalLoad) in prescribedLoads)
            {
                Node_v2[] nodes = FindBoundaryNodes(region, model, meshTolerance);
                double load = totalLoad / nodes.Length;
                foreach (Node_v2 node in nodes) model.Loads.Add(new Load_v2() { Node = node, DOF = dof, Amount = load });
            }

            return model;
        }

        public void PrescribeDisplacement(BoundaryRegion region, DOFType dof, double displacement)
            => prescribedDisplacements.Add((region, dof, displacement));

        /// <summary>
        /// </summary>
        /// <param name="load">Will be distributed evenly.</param>
        public void PrescribeDistributedLoad(BoundaryRegion region, DOFType dof, double load) 
            => prescribedLoads.Add((region, dof, load));

        private Node_v2[] FindBoundaryNodes(BoundaryRegion region, Model_v2 model, double tol)
        {
            double minX = 0.0, minY = 0.0, maxX = DomainLengthX, maxY = DomainLengthY; // for brevity

            IEnumerable<Node_v2> nodes;
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

            public bool Contains(CellConnectivity_v2 cell, double tol)
            {
                return cell.Vertices.All(node =>
                     (node.X >= minX - tol) && (node.X <= maxX + tol) && (node.Y >= minY - tol) && (node.Y <= maxY + tol));
            }
        }
    }
}
