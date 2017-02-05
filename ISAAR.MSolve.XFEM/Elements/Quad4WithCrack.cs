using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Enrichments;
using ISAAR.MSolve.XFEM.Enrichments.Jump;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Integration.GaussPoints;
using ISAAR.MSolve.XFEM.Integration.ShapeFunctions;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    class Crack2D : IEnrichmentFunction2D
    {
        private readonly ICurve2D curve;
        //private readonly PolynomialSmoothedHeaviside heaviside = new PolynomialSmoothedHeaviside(0.02);

        public Crack2D(ICurve2D curve)
        {
            this.curve = curve;
        }

        public double ValueAt(IPoint2D point)
        {
            double distance = curve.SignedDistanceOf(point);
            return (distance >= 0.0) ? 1.0 : 0.0;
        }

        public Tuple<double, double> DerivativesAt(IPoint2D point)
        {
            return new Tuple<double, double>(0.0, 0.0);
        }
    }

    class Quad4WithCrack
    {
        private const double DETERMINANT_TOLERANCE = 0.00000001;
        private const int STD_DOFS_COUNT = 8;
        private const int ENRICHED_DOFS_COUNT = 8;

        private readonly double halfLengthX;
        private readonly double halfLengthY;
        private readonly double centroidX;
        private readonly double centroidY;
        private readonly IReadOnlyList<Node2D> nodes;
        private readonly IReadOnlyList<GaussPoint2D> gaussPoints;
        private readonly IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materials;
        private readonly IEnrichmentFunction2D enrichment;

        private Dictionary<Node2D, double> nodalEnrichmentValues; // mutable

        public Quad4WithCrack(Node2D[] nodes, IFiniteElementMaterial2D material, ICurve2D discontinuity)
        {
            // TODO: Add checks here
            if ((nodes[0].X != nodes[3].X) || (nodes[1].X != nodes[2].X)
                || (nodes[0].Y != nodes[1].Y) || (nodes[2].Y != nodes[3].Y))
            {
                throw new ArgumentException("The local cartsian system is not aligned to the global cartesian system or "
                    + "may not be even a rectangle. Use an isoparametric quad 4 instead or check the order of the nodes");
            }

            this.halfLengthX = 0.5 * (nodes[1].X - nodes[0].X);
            this.centroidX = 0.5 * (nodes[0].X + nodes[1].X);
            this.halfLengthY = 0.5 * (nodes[3].Y - nodes[0].Y);
            this.centroidY = 0.5 * (nodes[0].Y + nodes[3].Y);

            this.nodes = new List<Node2D>(nodes);
            this.gaussPoints = FindGaussPoints(); // TODO: integration in an enriched element is much more complex

            var materialsDict = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in gaussPoints)
            {
                materialsDict[point] = material.Clone();
            }
            this.materials = materialsDict;

            this.enrichment = new Crack2D(discontinuity);
        }

        public SymmetricMatrix2D<double> BuildStdStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(8);
            foreach (var gaussPoint in gaussPoints) // TODO: remove the integration logic from the element class
            {
                // Calculate the necessary quantities for the integration
                Matrix2D<double> deformation = CalculateStdDeformationMatrix(gaussPoint);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Gauss integration at this point
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(thickness * gaussPoint.Weight);
                Debug.Assert(partial.Rows == 8);
                Debug.Assert(partial.Columns == 8);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(partial, stiffness);
            }
            return stiffness;
        }

        public void BuildEnrichedStiffnessMatrices(out Matrix2D<double> stiffnessStdEnriched,
            out SymmetricMatrix2D<double> stiffnessEnriched)
        {
            CalculateNodalEnrichmentValues();

            stiffnessStdEnriched = new Matrix2D<double>(STD_DOFS_COUNT, ENRICHED_DOFS_COUNT);
            stiffnessEnriched = new SymmetricMatrix2D<double>(ENRICHED_DOFS_COUNT);
            foreach (var gaussPoint in gaussPoints)
            {
                // Calculate the necessary quantities for the integration
                Matrix2D<double> Bstd = CalculateStdDeformationMatrix(gaussPoint);
                Matrix2D<double> Benr = CalculateEnrichedDeformationMatrix(gaussPoint);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Contributions of this gauss point to the element stiffness matrices
                Matrix2D<double> Kse = (Bstd.Transpose() * constitutive) * Benr;  // standard-enriched part
                Kse.Scale(thickness * gaussPoint.Weight);
                MatrixUtilities.AddPartialToTotalMatrix(Kse, stiffnessStdEnriched);

                Matrix2D<double> Kee = (Benr.Transpose() * constitutive) * Benr;  // enriched-enriched part
                Kee.Scale(thickness * gaussPoint.Weight);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(Kee, stiffnessEnriched);
            }
        }

        private void CalculateNodalEnrichmentValues()
        {
            nodalEnrichmentValues = new Dictionary<Node2D, double>();
            foreach (var node in nodes)
            {
                nodalEnrichmentValues[node] = enrichment.ValueAt(node);
            }
        }

        private IReadOnlyList<GaussPoint2D> FindGaussPoints()
        {
            var localGaussPoints = new List<GaussPoint2D>();
            foreach (var naturalGaussPoint in IntegrationRule2D.Order2x2.Points)
            {
                double localX = naturalGaussPoint.X * halfLengthX + centroidX;
                double localY = naturalGaussPoint.Y * halfLengthY + centroidY;
                double localWeight = naturalGaussPoint.Weight * halfLengthX * halfLengthY;
                localGaussPoints.Add(new GaussPoint2D(localX, localY, localWeight));
            }
            return localGaussPoints;
        }

        /// <summary>
        /// Calculate the deformation matrix B (3x8).
        /// B is a linear transformation FROM the nodal values of the displacement field TO the strain vector: 
        /// {e} = [B] * {d} => {u,x v,y u,y+v,x} = [B] * {u1 v1 u2 v2 u3 v3 u4 v4}
        /// </summary>
        /// <param name="gaussPoint"></param>
        /// <returns>A 3x8 matrix</returns>
        private Matrix2D<double> CalculateStdDeformationMatrix(GaussPoint2D gaussPoint)
        {
            var shapeFunctions = new Quad4ShapeFunctions(halfLengthX, halfLengthY);
            Tuple<double, double>[] shapeFunctionDerivatives =
                shapeFunctions.AllDerivativesAt(gaussPoint.X, gaussPoint.Y);

            var Bstd = new Matrix2D<double>(3, STD_DOFS_COUNT);
            for (int nodeIndex = 0; nodeIndex < 4; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                double Nx = shapeFunctionDerivatives[nodeIndex].Item1;
                double Ny = shapeFunctionDerivatives[nodeIndex].Item2;

                Bstd[0, col1] = Nx;
                Bstd[1, col2] = Ny;
                Bstd[2, col1] = Ny;
                Bstd[2, col2] = Nx;
            }
            return Bstd;
        }

        private Matrix2D<double> CalculateEnrichedDeformationMatrix(GaussPoint2D gaussPoint) 
        {
            // TODO: Evaluate the shape functions once and then pass the values, derivatives to the Bstd, Benr methods
            var shapeFunctions = new Quad4ShapeFunctions(halfLengthX, halfLengthY);
            double[] shapeFunctionValues = shapeFunctions.AllValuesAt(gaussPoint.X, gaussPoint.Y);
            var shapeFunctionDerivatives = shapeFunctions.AllDerivativesAt(gaussPoint.X, gaussPoint.Y);

            double F = enrichment.ValueAt(gaussPoint);
            var enrichmentDerivatives = enrichment.DerivativesAt(gaussPoint);
            double Fx = enrichmentDerivatives.Item1;
            double Fy = enrichmentDerivatives.Item2;

            var Benr = new Matrix2D<double>(3, 8); //TODO: Abstract the number of enriched dofs (8 here and the 2*node+1 afterwards)
            for (int nodeIndex = 0; nodeIndex < 4; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;

                double N = shapeFunctionValues[nodeIndex];
                double Nx = shapeFunctionDerivatives[nodeIndex].Item1;
                double Ny = shapeFunctionDerivatives[nodeIndex].Item2;
                double Fnode = nodalEnrichmentValues[nodes[nodeIndex]];

                Benr[0, col1] = Nx * (F - Fnode) + N * Fx;
                Benr[1, col2] = Ny * (F - Fnode) + N * Fy;
                Benr[2, col1] = Ny * (F - Fnode) + N * Fy;
                Benr[2, col2] = Nx * (F - Fnode) + N * Fx;
            }
            return Benr;
        }
    }
}
