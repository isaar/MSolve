using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    // TODO: Uses the same interpolation and nodes as the underlying std element! This must change
    // TODO: Enumerating artificial dofs may be needed to be done by this class. (e.g if structural FE introduce more artificial dofs than continuum FE)
    class XElement2D
    {
        public ContinuumElement2D StandardFiniteElement { get; }
        public IReadOnlyList<XNode2D> Nodes { get; } // The same as in stdElement.

        public static XElement2D CreateHomogeneous(ContinuumElement2D stdFiniteElement, 
            IFiniteElementMaterial2D commonMaterial)
        {
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> gpToMaterials =
                MaterialUtilities.AssignMaterialToIntegrationPoints(integration.GenerateIntegrationPoints(),
                commonMaterial);
            return new XElement2D(stdFiniteElement, gpToMaterials);
        }

        public static XElement2D CreateBimaterial(ContinuumElement2D stdFiniteElement,
            IFiniteElementMaterial2D materialLeft, IFiniteElementMaterial2D materialRight)
        {
            var integration = new SubgridIntegration2D(2, GaussQuadrature2D.Order2x2);
            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in integration.GenerateIntegrationPoints())
            {
                if (point.Xi < 0) gpToMaterials[point] = materialLeft.Clone();
                else gpToMaterials[point] = materialRight.Clone();
            }
            return new XElement2D(stdFiniteElement, gpToMaterials);
        }

        // Perhaps an element factory for the std element would be better than an instance,
        // in which I change stuff or use casts.
        private XElement2D(ContinuumElement2D stdFiniteElement,
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materialsOfGaussPoints)
        {
            this.StandardFiniteElement = stdFiniteElement;
            this.Nodes = (IReadOnlyList<XNode2D>)(stdFiniteElement.Nodes); // I am not too thrilled about casting especially when using covariance.
            this.StandardFiniteElement.MaterialsOfGaussPoints = materialsOfGaussPoints;
        }

        public SymmetricMatrix2D<double> BuildStdStiffnessMatrix() // This is not needed if I expose the standard FE
        {
            return StandardFiniteElement.BuildStiffnessMatrix();
        }

        #region Enriched stiffness matrices with Kse
        //public void BuildEnrichedStiffnessMatrices(out Matrix2D<double> stiffnessStdEnriched,
        //    out SymmetricMatrix2D<double> stiffnessEnriched)
        //{
        //    int artificialDofsCount = CountArtificialDofs();
        //    stiffnessStdEnriched = new Matrix2D<double>(StandardFiniteElement.DofsCount, artificialDofsCount);
        //    stiffnessEnriched = new SymmetricMatrix2D<double>(artificialDofsCount);
        //    foreach (var entry in StandardFiniteElement.MaterialsOfGaussPoints)
        //    {
        //        GaussPoint2D gaussPoint = entry.Key;
        //        IFiniteElementMaterial2D material = entry.Value;

        //        // Calculate the necessary quantities for the integration
        //        Matrix2D<double> constitutive = material.CalculateConstitutiveMatrix();
        //        EvaluatedInterpolation2D evaluatedInterpolation =
        //            StandardFiniteElement.Interpolation.EvaluateAt(Nodes, gaussPoint);
        //        Matrix2D<double> Bstd = StandardFiniteElement.CalculateDeformationMatrix(evaluatedInterpolation);
        //        Matrix2D<double> Benr = CalculateEnrichedDeformationMatrix(artificialDofsCount,
        //            gaussPoint, evaluatedInterpolation);

        //        // Contributions of this gauss point to the element stiffness matrices
        //        double dVolume = material.Thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight;
        //        Matrix2D<double> Kse = (Bstd.Transpose() * constitutive) * Benr;  // standard-enriched part
        //        Kse.Scale(dVolume); // Perhaps I should scale only the smallest matrix (constitutive) before the multiplications
        //        MatrixUtilities.AddPartialToTotalMatrix(Kse, stiffnessStdEnriched);

        //        Matrix2D<double> Kee = (Benr.Transpose() * constitutive) * Benr;  // enriched-enriched part
        //        Kee.Scale(dVolume);
        //        MatrixUtilities.AddPartialToSymmetricTotalMatrix(Kee, stiffnessEnriched);
        //    }
        //}
        #endregion

        public void BuildEnrichedStiffnessMatrices(out Matrix2D<double> stiffnessEnrichedStandard,
            out SymmetricMatrix2D<double> stiffnessEnriched)
        {
            int artificialDofsCount = CountArtificialDofs();
            stiffnessEnrichedStandard = new Matrix2D<double>(artificialDofsCount, StandardFiniteElement.DofsCount);
            stiffnessEnriched = new SymmetricMatrix2D<double>(artificialDofsCount);
            foreach (var entry in StandardFiniteElement.MaterialsOfGaussPoints)
            {
                GaussPoint2D gaussPoint = entry.Key;
                IFiniteElementMaterial2D material = entry.Value;

                // Calculate the necessary quantities for the integration
                Matrix2D<double> constitutive = material.CalculateConstitutiveMatrix();
                EvaluatedInterpolation2D evaluatedInterpolation =
                    StandardFiniteElement.Interpolation.EvaluateAt(Nodes, gaussPoint);
                Matrix2D<double> Bstd = StandardFiniteElement.CalculateDeformationMatrix(evaluatedInterpolation);
                Matrix2D<double> Benr = CalculateEnrichedDeformationMatrix(artificialDofsCount,
                    gaussPoint, evaluatedInterpolation);

                // Contributions of this gauss point to the element stiffness matrices. 
                // Kee = SUM(Benr^T * E * Benr * dV), Kes = SUM(Benr^T * E * Bstd * dV)
                double dVolume = material.Thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight;
                Matrix2D<double> transposeBenrTimesConstitutive = Benr.Transpose() * constitutive; // cache the result

                Matrix2D<double> Kes = transposeBenrTimesConstitutive * Bstd;  // enriched-standard part
                Kes.Scale(dVolume); // TODO: Scale only the smallest matrix (constitutive) before the multiplications. Probably requires a copy of the constitutive matrix.
                MatrixUtilities.AddPartialToTotalMatrix(Kes, stiffnessEnrichedStandard);

                Matrix2D<double> Kee = transposeBenrTimesConstitutive * Benr;  // enriched-enriched part
                Kee.Scale(dVolume);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(Kee, stiffnessEnriched);
            }
        }

        private Matrix2D<double> CalculateEnrichedDeformationMatrix(int artificialDofsCount,
            GaussPoint2D gaussPoint, EvaluatedInterpolation2D evaluatedInterpolation)
        {
            ICartesianPoint2D cartesianPoint = evaluatedInterpolation.TransformNaturalToCartesian(gaussPoint);
            var uniqueFunctions = new Dictionary<IEnrichmentFunction2D, EvaluatedFunction2D>();

            var deformationMatrix = new Matrix2D<double>(3, artificialDofsCount);
            int currentColumn = 0;
            foreach (XNode2D node in Nodes)
            {
                double N = evaluatedInterpolation.GetValueOf(node);
                var dNdx = evaluatedInterpolation.GetCartesianDerivativesOf(node);

                foreach (var enrichment in node.EnrichmentFunctions)
                {
                    IEnrichmentFunction2D enrichmentFunction = enrichment.Item1;
                    double nodalEnrichmentValue = enrichment.Item2;

                    // The enrichment function probably has been evaluated when processing a previous node. Avoid reevaluation.
                    EvaluatedFunction2D evaluatedEnrichment;
                    if (!(uniqueFunctions.TryGetValue(enrichmentFunction, out evaluatedEnrichment))) //Only search once
                    {
                        evaluatedEnrichment = enrichmentFunction.EvaluateAllAt(cartesianPoint);
                        uniqueFunctions[enrichmentFunction] = evaluatedEnrichment;
                    }

                    // For each node and with all derivatives w.r.t. cartesian coordinates, the enrichment derivatives 
                    // are: Bx = enrN,x = N,x(x,y) * [H(x,y) - H(node)] + N(x,y) * H,x(x,y), where H is the enrichment 
                    // function
                    double Bx = dNdx.Item1 * (evaluatedEnrichment.Value - nodalEnrichmentValue)
                        + N * evaluatedEnrichment.CartesianDerivatives.Item1;
                    double By = dNdx.Item2 * (evaluatedEnrichment.Value - nodalEnrichmentValue)
                        + N * evaluatedEnrichment.CartesianDerivatives.Item2;

                    // This depends on the convention: node major or enrichment major. The following is node major.
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

        public int CountArtificialDofs()
        {
            int count = 0;
            foreach (XNode2D node in Nodes) count += node.ArtificialDofsCount; // in all nodes or in enriched interpolation nodes?
            return count;
        }
    }
}
