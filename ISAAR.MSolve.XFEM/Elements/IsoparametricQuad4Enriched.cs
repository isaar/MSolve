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
    class IsoparametricQuad4Enriched
    {
        private const int NODES_COUNT = 4;
        private const int STD_DOFS_COUNT = 8;
        private const int ENRICHED_DOFS_COUNT = 8;

        private readonly IReadOnlyList<Node2D> nodes;
        private readonly IReadOnlyList<GaussPoint2D> gaussPoints;
        private readonly IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materials;
        private readonly IEnrichmentFunction2D enrichmentFunction;

        // TODO: use dictionaries for the nodal values.
        private double[] nodalEnrichmentValues; // mutable


        public IsoparametricQuad4Enriched(Node2D[] nodes, IFiniteElementMaterial2D material, 
            IEnrichmentFunction2D enrichmentFunction)
        {
            // TODO: Add checks here: order of nodes
            this.nodes = new List<Node2D>(nodes);
            //this.gaussPoints = IntegrationRule2D.Order10x10.Points; // TODO: remove the integration logic from the element class
            var integration = new SubgridIntegration(2, IntegrationRule2D.Order2x2);
            this.gaussPoints = integration.GeneratePoints();

            var materialsDict = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in gaussPoints)
            {
                materialsDict[point] = material.Clone();
            }
            this.materials = materialsDict;

            this.enrichmentFunction = enrichmentFunction;
        }

        public IsoparametricQuad4Enriched(Node2D[] nodes, IFiniteElementMaterial2D materialLeft, IFiniteElementMaterial2D materialRight,
            IEnrichmentFunction2D enrichmentFunction)
        {
            // TODO: Add checks here: order of nodes
            this.nodes = new List<Node2D>(nodes);
            //this.gaussPoints = IntegrationRule2D.Order10x10.Points; // TODO: remove the integration logic from the element class
            var integration = new SubgridIntegration(2, IntegrationRule2D.Order2x2);
            this.gaussPoints = integration.GeneratePoints();

            var materialsDict = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in gaussPoints)
            {
                if (point.X < 0) materialsDict[point] = materialLeft.Clone();
                else materialsDict[point] = materialRight.Clone();
            }
            this.materials = materialsDict;

            this.enrichmentFunction = enrichmentFunction;
        }

        public SymmetricMatrix2D<double> BuildStdStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(STD_DOFS_COUNT);
            foreach (var gaussPoint in gaussPoints)
            {
                // Calculate the necessary quantities for the integration
                ShapeFunctionDerivatives2D shapeFunctionDerivatives =
                    IsoparametricQuad4ShapeFunctions.AllDerivativesAt(gaussPoint.X, gaussPoint.Y);
                XJacobian2D jacobian = new XJacobian2D(nodes, shapeFunctionDerivatives);

                Matrix2D<double> deformation = CalculateStdDeformationMatrix(shapeFunctionDerivatives, jacobian);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Gauss integration at this point
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(thickness * jacobian.Determinant * gaussPoint.Weight);
                Debug.Assert(partial.Rows == STD_DOFS_COUNT);
                Debug.Assert(partial.Columns == STD_DOFS_COUNT);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(partial, stiffness);
            }
            return stiffness;
        }

        public void BuildEnrichedStiffnessMatrices(out Matrix2D<double> stiffnessStdEnriched,
            out SymmetricMatrix2D<double> stiffnessEnriched)
        {
            CalculateNodalEnrichments();

            stiffnessStdEnriched = new Matrix2D<double>(STD_DOFS_COUNT, ENRICHED_DOFS_COUNT);
            stiffnessEnriched = new SymmetricMatrix2D<double>(ENRICHED_DOFS_COUNT);
            foreach (var gaussPoint in gaussPoints)
            {
                // Calculate the necessary quantities for the integration
                ShapeFunctionDerivatives2D shapeFunctionDerivatives =
                    IsoparametricQuad4ShapeFunctions.AllDerivativesAt(gaussPoint.X, gaussPoint.Y);
                XJacobian2D jacobian = new XJacobian2D(nodes, shapeFunctionDerivatives);

                Matrix2D<double> Bstd = CalculateStdDeformationMatrix(shapeFunctionDerivatives, jacobian);
                Matrix2D<double> Benr = CalculateEnrichedDeformationMatrix(gaussPoint, shapeFunctionDerivatives, jacobian);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Contributions of this gauss point to the element stiffness matrices
                double scalar = thickness * jacobian.Determinant * gaussPoint.Weight;
                Matrix2D<double> Kse = (Bstd.Transpose() * constitutive) * Benr;  // standard-enriched part
                Kse.Scale(scalar);
                MatrixUtilities.AddPartialToTotalMatrix(Kse, stiffnessStdEnriched);

                Matrix2D<double> Kee = (Benr.Transpose() * constitutive) * Benr;  // enriched-enriched part
                Kee.Scale(scalar);
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

        private Matrix2D<double> CalculateStdDeformationMatrix(ShapeFunctionDerivatives2D shapeFunctionDerivatives,
            XJacobian2D jacobian)
        {
            //Calculate B2 deformation matrix. Dimensions = 4x8.
            //B2 is a linear transformation FROM the nodal values of the displacement field TO the the derivatives of
            //the displacement field in respect to the natural axes: {dU/dXi} = [B2] * {d} => 
            //{u,xi u,eta v,xi, v,eta} = [B2] * {u1 v1 u2 v2 u3 v3 u4 v4}
            var B2 = new Matrix2D<double>(4, STD_DOFS_COUNT);
            for (int nodeIndex = 0; nodeIndex < 4; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                double Nxi = shapeFunctionDerivatives.XiDerivativeOfNode(nodeIndex);
                double Neta = shapeFunctionDerivatives.EtaDerivativeOfNode(nodeIndex);

                B2[0, col1] = Nxi;
                B2[1, col1] = Neta;
                B2[2, col2] = Nxi;
                B2[3, col2] = Neta;
            }

            // The deformation matrix can be calculated as [B] = [B1] * [B2]
            return jacobian.CalculateB1DeformationMatrix() * B2;
        }

        private Matrix2D<double> CalculateEnrichedDeformationMatrix(GaussPoint2D gaussPoint,
            ShapeFunctionDerivatives2D shapeFunctionDerivatives, XJacobian2D jacobian)
        {
            // TODO: Evaluate the shape functions only ONCE and then pass the values, derivatives to the Bstd, Benr methods
            double[] shapeFunctionValues = IsoparametricQuad4ShapeFunctions.AllValuesAt(gaussPoint.X, gaussPoint.Y);

            double x = InterpolationUtilities.InterpolateNodalValuesToPoint(NodalX, shapeFunctionValues);
            double y = InterpolationUtilities.InterpolateNodalValuesToPoint(NodalY, shapeFunctionValues);
            Point2D cartesianPoint = new Point2D(x, y);

            // Calculate the enrichment value and derivatives at this Gauss point
            double enrichmentValue = enrichmentFunction.ValueAt(cartesianPoint);
            var enrichmentDerivativesCartesian = enrichmentFunction.DerivativesAt(cartesianPoint);
            double enrichmentDerivativeXi = enrichmentDerivativesCartesian.Item1 * jacobian.dXdXi + 
                enrichmentDerivativesCartesian.Item2 * jacobian.dYdXi;
            double enrichmentDerivativeEta = enrichmentDerivativesCartesian.Item1 * jacobian.dXdEta +
                enrichmentDerivativesCartesian.Item2 * jacobian.dYdEta;

            var B2enr = new Matrix2D<double>(4, 8); //TODO: Abstract the number of enriched dofs (8 here and the 2*node+1 afterwards)
            for (int nodeIndex = 0; nodeIndex < 4; ++nodeIndex)
            {
                double N = shapeFunctionValues[nodeIndex];
                double Nxi = shapeFunctionDerivatives.XiDerivativeOfNode(nodeIndex);
                double Neta = shapeFunctionDerivatives.EtaDerivativeOfNode(nodeIndex);
                double nodalEnrichment = nodalEnrichmentValues[nodeIndex];

                double enrN_xi = Nxi * (enrichmentValue - nodalEnrichment) + N * enrichmentDerivativeXi;
                double enrN_eta = Neta * (enrichmentValue - nodalEnrichment) + N * enrichmentDerivativeEta;

                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                B2enr[0, col1] = enrN_xi;
                B2enr[1, col1] = enrN_eta;
                B2enr[2, col2] = enrN_xi;
                B2enr[3, col2] = enrN_eta;
            }
            return jacobian.CalculateB1DeformationMatrix() * B2enr;
        }

        private double[] NodalX
        {
            get
            {
                double[] x = new double[NODES_COUNT];
                for (int nodeId = 0; nodeId < NODES_COUNT; ++nodeId) x[nodeId] = nodes[nodeId].X;
                return x;
            }
        }

        private double[] NodalY
        {
            get
            {
                double[] y = new double[NODES_COUNT];
                for (int nodeId = 0; nodeId < NODES_COUNT; ++nodeId) y[nodeId] = nodes[nodeId].Y;
                return y;
            }
        }
    }
}
