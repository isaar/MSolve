using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
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
    // TODO: It still uses the same nodes for the std and enriched part. 
    class XIsoparametricQuad4
    {
        //private const int NODES_COUNT = 4;
        //private const int ENRICHMENT_FUNCTIONS_COUNT = 4;
        //private const int STD_DOFS_COUNT = 8;
        //private const int ENRICHED_DOFS_COUNT = 32;
        //private const int ENRICHED_DOFS_PER_NODE = 8;

        // Copy of the nodes list in the std FE? Or only the enrichment interpolation nodes? Should I use the list stored in the enrichment interpolation?
        public IReadOnlyList<XNode2D> Nodes { get; } 
        private readonly IsoparametricQuad4 stdFiniteElement;

        // I could store the materials and gauss points here (like the nodes), instead of pulling them from the std FE.
        //public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> MaterialsOfGaussPoints { get; } 

        private readonly IsoparametricInterpolation2D enrichmentInterpolation;
        

        public static XIsoparametricQuad4 CreateHomogeneous(XNode2D[] nodes, IFiniteElementMaterial2D material)
        {
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in integration.GenerateIntegrationPoints())
            {
                gpToMaterials[point] = material.Clone();
            }
            return new XIsoparametricQuad4(nodes, gpToMaterials);
        }

        public static XIsoparametricQuad4 CreateBimaterial(XNode2D[] nodes, 
            IFiniteElementMaterial2D materialLeft, IFiniteElementMaterial2D materialRight)
        {
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in integration.GenerateIntegrationPoints())
            {
                if (point.X < 0) gpToMaterials[point] = materialLeft.Clone();
                else gpToMaterials[point] = materialRight.Clone();
            }
            return new XIsoparametricQuad4(nodes, gpToMaterials);
        }

        private XIsoparametricQuad4(XNode2D[] nodes, 
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materialsOfGaussPoints)
        {
            var nodesCopy = new XNode2D[nodes.Length];
            nodes.CopyTo(nodesCopy, 0);
            this.Nodes = nodesCopy;

            // Checking the nodes and gauss points is done by the standard Finite Element
            // As it is, the same nodes are used for both the std FE and the enriched FE.
            this.stdFiniteElement = new IsoparametricQuad4(this.Nodes, materialsOfGaussPoints);
            this.enrichmentInterpolation = new IsoparametricInterpolation2D(this.Nodes, NaturalShapeFunctions2D.Quad4);
        }

        public SymmetricMatrix2D<double> BuildStdStiffnessMatrix()
        {
            return stdFiniteElement.BuildStiffnessMatrix();
        }

        public void BuildEnrichedStiffnessMatrices(out Matrix2D<double> stiffnessStdEnriched,
            out SymmetricMatrix2D<double> stiffnessEnriched)
        {
            int artificialDofsCount = CountArtificialDofs();
            stiffnessStdEnriched = new Matrix2D<double>(stdFiniteElement.DOFS_COUNT, artificialDofsCount);
            stiffnessEnriched = new SymmetricMatrix2D<double>(artificialDofsCount);
            foreach (var entry in stdFiniteElement.MaterialsOfGaussPoints)
            {
                GaussPoint2D gaussPoint = entry.Key;
                IFiniteElementMaterial2D material = entry.Value;

                // Calculate the necessary quantities for the integration
                Matrix2D<double> constitutive = material.CalculateConstitutiveMatrix();
                EvaluatedInterpolation2D evaluatedInterpolation = enrichmentInterpolation.EvaluateAt(gaussPoint);
                Matrix2D<double> Bstd = stdFiniteElement.CalculateDeformationMatrix(evaluatedInterpolation);
                Matrix2D<double> Benr = CalculateEnrichedDeformationMatrix(artificialDofsCount, 
                    gaussPoint, evaluatedInterpolation);

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

        private Matrix2D<double> CalculateEnrichedDeformationMatrix(int artificialDofsCount, 
            GaussPoint2D gaussPoint, EvaluatedInterpolation2D evaluatedInterpolation)
        {
            IPoint2D cartesianPoint = evaluatedInterpolation.TransformNaturalToCartesian(gaussPoint);

            // Calculate the enrichment values and derivatives at this Gauss point. 
            // TODO: If the enrichment functions are separate for each node, their common parts must be calculated only
            // once. This is probably the step to do it.
            //double[] enrichmentValues = new double[ENRICHMENT_FUNCTIONS_COUNT];
            //var enrichmentDerivatives = new Tuple<double, double>[ENRICHMENT_FUNCTIONS_COUNT];
            //for (int enrichmentIdx = 0; enrichmentIdx < ENRICHMENT_FUNCTIONS_COUNT; ++enrichmentIdx)
            //{
            //    enrichmentValues[enrichmentIdx] = enrichmentFunctions[enrichmentIdx].ValueAt(cartesianPoint);
            //    enrichmentDerivatives[enrichmentIdx] = enrichmentFunctions[enrichmentIdx].DerivativesAt(cartesianPoint);
            //}

            // Build the deformation matrix per node.
            var deformationMatrix = new Matrix2D<double>(3, artificialDofsCount);
            int currentColumn = 0;
            foreach (XNode2D node in Nodes)
            {
                double N = evaluatedInterpolation.GetValueOf(node);
                var dNdx = evaluatedInterpolation.GetCartesianDerivativesOf(node);

                foreach (var enrichment in node.EnrichmentFunctions)
                {
                    // Denote the enrichment function as H
                    double H = enrichment.Item1.ValueAt(cartesianPoint); // TODO: Calculate the values and derivatives of all unique enrichment functions only once, instead of once per associated node.
                    Tuple<double, double> dHdx = enrichment.Item1.DerivativesAt(cartesianPoint);
                    double nodalH = enrichment.Item2;

                    // For each node and with all derivatives w.r.t. cartesian coordinates, the enrichment derivatives are:
                    // Bx = enrN,x = N,x(x,y) * [H(x,y) - H(node)] + N(x,y) * H,x(x,y). Similarly for By
                    double Bx = dNdx.Item1 * (H - nodalH) + N * dHdx.Item1;
                    double By = dNdx.Item2 * (H - nodalH) + N * dHdx.Item2;

                    // This depends on the convention: node major or enrichment major. The following is node major
                    int col1 = currentColumn++;
                    int col2 = currentColumn++;

                    deformationMatrix[0, col1] = Bx;
                    deformationMatrix[1, col2] = By;
                    deformationMatrix[2, col1] = By;
                    deformationMatrix[2, col2] = Bx;
                }
            }
            Debug.Assert(currentColumn == artificialDofsCount);
            return deformationMatrix;
        }

        private int CountArtificialDofs()
        {
            int count = 0;
            foreach (XNode2D node in Nodes) count += node.ArtificialDofsCount; // in all nodes or in enriched interpolation nodes?
            return count;
        }
    }
}
