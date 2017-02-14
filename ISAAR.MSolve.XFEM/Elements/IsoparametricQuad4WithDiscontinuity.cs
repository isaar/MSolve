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

        private readonly IReadOnlyList<Node2D> nodes; // Copy of the nodes list in the std FE
        private readonly IsoparametricQuad4 stdFiniteElement;
        
        private readonly IsoparametricInterpolation2D enrichedInterpolation;
        private readonly IEnrichmentFunction2D enrichmentFunction;

        // I could store the materials and gauss points here (like the nodes), instead of pulling them from the std FE.
        //public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> MaterialsOfGaussPoints { get; } 

        // TODO: use dictionaries for the nodal values.
        private double[] nodalEnrichmentValues; // mutable

        public static IsoparametricQuad4WithDiscontinuity CreateHomogeneous(Node2D[] nodes, 
            IEnrichmentFunction2D enrichmentFunction, IFiniteElementMaterial2D material)
        {
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in integration.GenerateIntegrationPoints())
            {
                gpToMaterials[point] = material.Clone();
            }
            return new IsoparametricQuad4WithDiscontinuity(nodes, enrichmentFunction, gpToMaterials);
        }

        public static IsoparametricQuad4WithDiscontinuity CreateBimaterial(Node2D[] nodes, IEnrichmentFunction2D enrichmentFunction,
            IFiniteElementMaterial2D materialLeft, IFiniteElementMaterial2D materialRight)
        {
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in integration.GenerateIntegrationPoints())
            {
                if (point.X < 0) gpToMaterials[point] = materialLeft.Clone();
                else gpToMaterials[point] = materialRight.Clone();
            }
            return new IsoparametricQuad4WithDiscontinuity(nodes, enrichmentFunction, gpToMaterials);
        }

        private IsoparametricQuad4WithDiscontinuity(Node2D[] nodes, IEnrichmentFunction2D enrichmentFunction,
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materialsOfGaussPoints)
        {
            var nodesCopy = new Node2D[nodes.Length];
            nodes.CopyTo(nodesCopy, 0);
            this.nodes = nodesCopy;
            // Checking the nodes and gauss points is done by the standard Finite Element
            this.stdFiniteElement = new IsoparametricQuad4(this.nodes, materialsOfGaussPoints); 
            this.enrichedInterpolation = new IsoparametricInterpolation2D(this.nodes, NaturalShapeFunctions2D.Quad4);
            this.enrichmentFunction = enrichmentFunction;
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
                EvaluatedInterpolation2D evaluatedInterpolation = this.enrichedInterpolation.EvaluateAt(gaussPoint);
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
            nodalEnrichmentValues = new double[NODES_COUNT];
            for (int nodeIndex = 0; nodeIndex < NODES_COUNT; ++nodeIndex)
            {
                nodalEnrichmentValues[nodeIndex] = enrichmentFunction.EvalueAt(nodes[nodeIndex]);
            }
        }

        private Matrix2D<double> CalculateEnrichedDeformationMatrix(GaussPoint2D gaussPoint,
            EvaluatedInterpolation2D evaluatedInterpolation)
        {
            IPoint2D cartesianPoint = evaluatedInterpolation.TransformNaturalToCartesian(gaussPoint);

            // Calculate the enrichment value and derivatives at this Gauss point. Denote the enrichment function as H
            double H = enrichmentFunction.EvalueAt(cartesianPoint);
            Tuple<double, double> dHdx = enrichmentFunction.EvaluateDerivativesAt(cartesianPoint);
            
            // Build the deformation matrix per node
            var deformationMatrix = new Matrix2D<double>(3, 8); //TODO: Abstract the number of enriched dofs (8 here and the 2*node+1 afterwards)
            for (int nodeIdx = 0; nodeIdx < 4; ++nodeIdx)
            {
                Node2D node = nodes[nodeIdx];

                double N = evaluatedInterpolation.GetValueOf(node);
                Tuple<double, double> dNdx = evaluatedInterpolation.GetCartesianDerivativesOf(node);
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
