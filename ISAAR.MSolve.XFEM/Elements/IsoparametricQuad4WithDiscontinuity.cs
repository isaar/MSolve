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
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Interpolation.ShapeFunctions;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    class IsoparametricQuad4WithDiscontinuity
    {
        private const int NODES_COUNT = 4;
        private const int STD_DOFS_COUNT = 8;
        private const int ENRICHED_DOFS_COUNT = 8;

        private readonly IReadOnlyList<Node2D> nodes;
        private readonly IReadOnlyList<GaussPoint2D> gaussPoints;
        private readonly IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materials;
        private readonly IsoparametricInterpolation2D interpolation;
        private readonly IEnrichmentFunction2D enrichmentFunction;

        // TODO: use dictionaries for the nodal values.
        private double[] nodalEnrichmentValues; // mutable


        public IsoparametricQuad4WithDiscontinuity(Node2D[] nodes, IFiniteElementMaterial2D material, 
            IEnrichmentFunction2D enrichmentFunction)
        {
            // TODO: Add checks here: order of nodes
            this.nodes = new List<Node2D>(nodes);
            //this.gaussPoints = IntegrationRule2D.Order10x10.Points; // TODO: remove the integration logic from the element class
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            this.gaussPoints = integration.GenerateIntegrationPoints();

            var materialsDict = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in gaussPoints)
            {
                materialsDict[point] = material.Clone();
            }
            this.materials = materialsDict;

            this.interpolation = new IsoparametricInterpolation2D(this.nodes, NaturalShapeFunctions2D.Quad4);

            this.enrichmentFunction = enrichmentFunction;
        }

        public IsoparametricQuad4WithDiscontinuity(Node2D[] nodes, IFiniteElementMaterial2D materialLeft, IFiniteElementMaterial2D materialRight,
            IEnrichmentFunction2D enrichmentFunction)
        {
            // TODO: Add checks here: order of nodes
            this.nodes = new List<Node2D>(nodes);
            //this.gaussPoints = IntegrationRule2D.Order10x10.Points; // TODO: remove the integration logic from the element class
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            this.gaussPoints = integration.GenerateIntegrationPoints();

            var materialsDict = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in gaussPoints)
            {
                if (point.X < 0) materialsDict[point] = materialLeft.Clone();
                else materialsDict[point] = materialRight.Clone();
            }
            this.materials = materialsDict;

            this.interpolation = new IsoparametricInterpolation2D(this.nodes, NaturalShapeFunctions2D.Quad4);

            this.enrichmentFunction = enrichmentFunction;
        }

        public SymmetricMatrix2D<double> BuildStdStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(STD_DOFS_COUNT);
            foreach (var gaussPoint in gaussPoints) // These gauss points could only be the std 2x2 ones. How would materials work then?
            {
                // Calculate the necessary quantities for the integration
                Interpolation.ShapeFunctionDerivatives2D interpolationDerivatives =
                    this.interpolation.EvaluateDerivativesAt(gaussPoint.X, gaussPoint.Y);
                Matrix2D<double> deformation = CalculateStdDeformationMatrix(interpolationDerivatives);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(thickness * interpolationDerivatives.Jacobian.Determinant * gaussPoint.Weight);
                Debug.Assert(partial.Rows == STD_DOFS_COUNT);
                Debug.Assert(partial.Columns == STD_DOFS_COUNT);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(partial, stiffness);
            }
            return stiffness;
        }

        public void BuildEnrichedStiffnessMatrices(out Matrix2D<double> stiffnessStdEnriched,
            out SymmetricMatrix2D<double> stiffnessEnriched)
        {
            CalculateNodalEnrichments(); // Only needs to be done once for all gauss points

            stiffnessStdEnriched = new Matrix2D<double>(STD_DOFS_COUNT, ENRICHED_DOFS_COUNT);
            stiffnessEnriched = new SymmetricMatrix2D<double>(ENRICHED_DOFS_COUNT);
            foreach (var gaussPoint in gaussPoints)
            {
                // Calculate the necessary quantities for the integration
                Interpolation.ShapeFunctionDerivatives2D interpolationDerivatives =
                    this.interpolation.EvaluateDerivativesAt(gaussPoint.X, gaussPoint.Y);
                Matrix2D<double> Bstd = CalculateStdDeformationMatrix(interpolationDerivatives);
                Matrix2D<double> Benr = CalculateEnrichedDeformationMatrix(gaussPoint, interpolationDerivatives);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Contributions of this gauss point to the element stiffness matrices
                double dVolume = thickness * interpolationDerivatives.Jacobian.Determinant * gaussPoint.Weight;
                Matrix2D<double> Kse = (Bstd.Transpose() * constitutive) * Benr;  // standard-enriched part
                Kse.Scale(dVolume);
                MatrixUtilities.AddPartialToTotalMatrix(Kse, stiffnessStdEnriched);

                Matrix2D<double> Kee = (Benr.Transpose() * constitutive) * Benr;  // enriched-enriched part
                Kee.Scale(dVolume);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(Kee, stiffnessEnriched);
            }
        }

        private void CalculateNodalEnrichments()
        {
            nodalEnrichmentValues = new double[NODES_COUNT];
            for (int nodeIndex = 0; nodeIndex < NODES_COUNT; ++nodeIndex)
            {
                nodalEnrichmentValues[nodeIndex] = enrichmentFunction.ValueAt(nodes[nodeIndex]);
            }
        }

        /// <summary>
        /// Calculates the deformation matrix B. Dimensions = 3x8.
        /// B is a linear transformation FROM the nodal values of the displacement field TO the the derivatives of
        /// the displacement field in respect to the cartesian axes (i.e. the stresses): {dU/dX} = [B] * {d} => 
        /// {u,x v,y u,y, v,x} = [... Bk ...] * {u1 v1 u2 v2 u3 v3 u4 v4}, where k = 1, ... nodesCount is a node and
        /// Bk = [dNk/dx 0; 0 dNk/dY; dNk/dy dNk/dx] (3x2)
        /// </summary>
        /// <param name="shapeFunctionDerivatives">The shape function derivatives calculated at a specific 
        ///     integration point</param>
        /// <returns></returns>
        private Matrix2D<double> CalculateStdDeformationMatrix(ShapeFunctionDerivatives2D shapeFunctionDerivatives)
        {
            var deformationMatrix = new Matrix2D<double>(3, STD_DOFS_COUNT);
            for (int nodeIndex = 0; nodeIndex < 4; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                Tuple<double, double> dNdx = shapeFunctionDerivatives.CartesianDerivativesOfNode(nodeIndex);

                deformationMatrix[0, col1] = dNdx.Item1;
                deformationMatrix[1, col2] = dNdx.Item2;
                deformationMatrix[2, col1] = dNdx.Item2;
                deformationMatrix[2, col2] = dNdx.Item1;
            }
            return deformationMatrix;
        }

        private Matrix2D<double> CalculateEnrichedDeformationMatrix(GaussPoint2D gaussPoint,
            ShapeFunctionDerivatives2D shapeFunctionDerivatives)
        {
            ShapeFunctionValues2D shapeFunctionValues = interpolation.EvaluateAt(gaussPoint.X, gaussPoint.Y);
            IPoint2D cartesianPoint = shapeFunctionValues.TransformNaturalToCartesian(gaussPoint);

            // Calculate the enrichment value and derivatives at this Gauss point. Denote the enrichment function as H
            double H = enrichmentFunction.ValueAt(cartesianPoint);
            Tuple<double, double> dHdx = enrichmentFunction.DerivativesAt(cartesianPoint);
            
            // Build the deformation matrix per node
            var deformationMatrix = new Matrix2D<double>(3, 8); //TODO: Abstract the number of enriched dofs (8 here and the 2*node+1 afterwards)
            for (int nodeIdx = 0; nodeIdx < 4; ++nodeIdx)
            {
                double N = shapeFunctionValues[nodeIdx];
                Tuple<double, double> dNdx = shapeFunctionDerivatives.CartesianDerivativesOfNode(nodeIdx);
                double nodalH = nodalEnrichmentValues[nodeIdx];

                // For each node and with all derivatives w.r.t. cartesian coordinates, the enrichment derivatives are:
                // Bx = enrN,x = N,x(x,y) * [enrichment(x,y) - enrichment(node)] + N(x,y) * enrichment,x(x,y)
                double Bx = dNdx.Item1 * (H - nodalH) + N * dHdx.Item1;
                double By = dNdx.Item2 * (H - nodalH) + N * dHdx.Item2;

                int col1 = 2 * nodeIdx;
                int col2 = 2 * nodeIdx + 1;
                deformationMatrix[0, col1] = Bx;
                deformationMatrix[1, col2] = By;
                deformationMatrix[2, col1] = By;
                deformationMatrix[2, col2] = Bx;
            }
            return deformationMatrix;
        }
    }
}
