using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    // dofs = {u d1 d2 d3 d4}, where: 
    // u =  {u1 v1 u2 v2 u3 v3 u4 v4} (8x1) are the standard displacement dofs
    // d1...d4 are all the artificial dofs for the nodes 1...4 respectively. In more detail for the node k: 
    // dk = {b1xk b1yk b2xk b2yk b3xk b3yk b4xk b4yk} (8x1), 
    // where bixk is the x dof corresponding to the ith tip enrichment function that has the support of node k
    // Total enriched dofs = 8 artificial dofs/node * 4 nodes = 32
    class IsoparametricQuad4WithTip
    {
        private const int NODES_COUNT = 4;
        private const int ENRICHMENT_FUNCTIONS_COUNT = 4;
        private const int STD_DOFS_COUNT = 8;
        private const int ENRICHED_DOFS_COUNT = 32;
        private const int ENRICHED_DOFS_PER_NODE = 8;

        private readonly IReadOnlyList<Node2D> nodes; // Copy of the nodes list in the std FE
        private readonly IsoparametricQuad4_OLD stdFiniteElement;

        // I could store the materials and gauss points here (like the nodes), instead of pulling them from the std FE.
        //public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> MaterialsOfGaussPoints { get; } 

        private readonly IsoparametricInterpolation2D enrichedInterpolation;
        // TODO: use dictionaries for the nodal values and functions.
        private readonly IEnrichmentFunction2D[] enrichmentFunctions;
        private double[,] nodalEnrichmentValues; // mutable

        public static IsoparametricQuad4WithTip CreateHomogeneous(Node2D[] nodes,
            IEnrichmentFunction2D[] enrichmentFunctions, IFiniteElementMaterial2D material)
        {
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in integration.GenerateIntegrationPoints())
            {
                gpToMaterials[point] = material.Clone();
            }
            return new IsoparametricQuad4WithTip(nodes, enrichmentFunctions, gpToMaterials);
        }

        private IsoparametricQuad4WithTip(Node2D[] nodes, IEnrichmentFunction2D[] enrichmentFunctions,
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materialsOfGaussPoints)
        {
            var nodesCopy = new Node2D[nodes.Length];
            nodes.CopyTo(nodesCopy, 0);
            this.nodes = nodesCopy;
            // Checking the nodes and gauss points is done by the standard Finite Element
            this.stdFiniteElement = new IsoparametricQuad4_OLD(this.nodes, materialsOfGaussPoints);
            this.enrichedInterpolation = IsoparametricInterpolation2D.Quad4;
            this.enrichmentFunctions = enrichmentFunctions;
        }

        public SymmetricMatrix2D<double> BuildStdStiffnessMatrix()
        {
            return stdFiniteElement.BuildStiffnessMatrix();
        }

        public void BuildEnrichedStiffnessMatrices(out Matrix2D<double> stiffnessStdEnriched,
            out SymmetricMatrix2D<double> stiffnessEnriched)
        {
            CalculateNodalEnrichments(); // Only needs to be done once for all gauss points

            stiffnessStdEnriched = new Matrix2D<double>(STD_DOFS_COUNT, ENRICHED_DOFS_COUNT);
            stiffnessEnriched = new SymmetricMatrix2D<double>(ENRICHED_DOFS_COUNT);
            foreach (var entry in stdFiniteElement.MaterialsOfGaussPoints)
            {
                GaussPoint2D gaussPoint = entry.Key;
                IFiniteElementMaterial2D material = entry.Value;

                // Calculate the necessary quantities for the integration
                EvaluatedInterpolation2D evaluatedInterpolation = enrichedInterpolation.EvaluateAt(nodes, gaussPoint);
                Matrix2D<double> Bstd = stdFiniteElement.CalculateDeformationMatrix(evaluatedInterpolation);
                Matrix2D<double> Benr = CalculateEnrichedDeformationMatrix(gaussPoint, evaluatedInterpolation);
                Matrix2D<double> constitutive = material.CalculateConstitutiveMatrix();

                // Contributions of this gauss point to the element stiffness matrices
                double dVolume = material.Thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight;
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
            nodalEnrichmentValues = new double[NODES_COUNT, ENRICHMENT_FUNCTIONS_COUNT];
            for (int nodeIdx = 0; nodeIdx < NODES_COUNT; ++nodeIdx)
            {
                for (int enrichmentIdx = 0; enrichmentIdx < ENRICHMENT_FUNCTIONS_COUNT; ++enrichmentIdx)
                {
                    nodalEnrichmentValues[nodeIdx, enrichmentIdx] = 
                        enrichmentFunctions[enrichmentIdx].EvalueAt(nodes[nodeIdx]);
                }
            }
        }

        private Matrix2D<double> CalculateEnrichedDeformationMatrix(GaussPoint2D gaussPoint,
            EvaluatedInterpolation2D evaluatedInterpolation)
        {
            IPoint2D cartesianPoint = evaluatedInterpolation.TransformNaturalToCartesian(gaussPoint);

            // Calculate the enrichment values and derivatives at this Gauss point. 
            // TODO: If the enrichment functions are separate for each node, their common parts must be calculated only
            // once. This is probably the step to do it.
            double[] enrichmentValues = new double[ENRICHMENT_FUNCTIONS_COUNT];
            var enrichmentDerivatives = new Tuple<double, double>[ENRICHMENT_FUNCTIONS_COUNT];
            for (int enrichmentIdx = 0; enrichmentIdx < ENRICHMENT_FUNCTIONS_COUNT; ++enrichmentIdx)
            {
                enrichmentValues[enrichmentIdx] = enrichmentFunctions[enrichmentIdx].EvalueAt(cartesianPoint);
                enrichmentDerivatives[enrichmentIdx] = enrichmentFunctions[enrichmentIdx].EvaluateDerivativesAt(cartesianPoint);
            }

            // Build the deformation matrix per node.
            var deformationMatrix = new Matrix2D<double>(3, ENRICHED_DOFS_COUNT); //TODO: Abstract the number of enriched dofs (8 here and the 2*node+1 afterwards)
            for (int nodeIdx = 0; nodeIdx < 4; ++nodeIdx)
            {
                Node2D node = nodes[nodeIdx];

                double N = evaluatedInterpolation.GetValueOf(node);
                Tuple<double, double> dNdx = evaluatedInterpolation.GetCartesianDerivativesOf(node);

                for (int enrichmentIdx = 0; enrichmentIdx < ENRICHMENT_FUNCTIONS_COUNT; ++enrichmentIdx)
                {
                    // Denote the enrichment function as H
                    double nodalH = nodalEnrichmentValues[nodeIdx, enrichmentIdx];
                    double H = enrichmentValues[enrichmentIdx];
                    Tuple<double, double> dHdx = enrichmentDerivatives[enrichmentIdx];

                    // For each node and with all derivatives w.r.t. cartesian coordinates, the enrichment derivatives are:
                    // Bx = enrN,x = N,x(x,y) * [H(x,y) - H(node)] + N(x,y) * H,x(x,y). Similarly for By
                    double Bx = dNdx.Item1 * (H - nodalH) + N * dHdx.Item1;
                    double By = dNdx.Item2 * (H - nodalH) + N * dHdx.Item2;

                    // This depends on the convention: node major or enrichment major. The following is node major
                    int col1 = ENRICHED_DOFS_PER_NODE * nodeIdx + 2 * enrichmentIdx;
                    int col2 = ENRICHED_DOFS_PER_NODE * nodeIdx + 2 * enrichmentIdx + 1;

                    deformationMatrix[0, col1] = Bx;
                    deformationMatrix[1, col2] = By;
                    deformationMatrix[2, col1] = By;
                    deformationMatrix[2, col2] = Bx;
                }
            }
            return deformationMatrix;
        }
    }
}
