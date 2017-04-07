using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Elements
{
    /// <summary>
    /// TODO: Uses the same interpolation and nodes as the underlying std element! This must change
    /// TODO: Enumerating artificial dofs may be needed to be done by this class. (e.g if structural FE introduce more 
    ///     artificial dofs than continuum FE)
    /// TODO: The calculation of Kss uses the same gauss points as the calculation of Kes, Kee. 
    ///     Pros: only need to track one set of Gauss points, which simplifies non linear analysis. 
    ///     Cons: calculating Kss with the Gauss points of an enriched element is much more expensive
    /// </summary>
    class XContinuumElement2D: IMeshFace
    {
        protected readonly IsoparametricElementType2D elementType;

        public IReadOnlyList<ICartesianPoint2D> Vertices { get { return Nodes; } }

        /// <summary>
        /// All nodes are enriched for now.
        /// </summary>
        public IReadOnlyList<XNode2D> Nodes { get; }

        /// <summary>
        /// Common interpolation for standard and enriched nodes.
        /// </summary>
        public IsoparametricInterpolation2D Interpolation { get { return elementType.Interpolation; } }

        public IStandardQuadrature2D StandardQuadrature { get { return elementType.StandardQuadrature; } }
        public IIntegrationStrategy2D<XContinuumElement2D> IntegrationStrategy { get; }

        /// <summary>
        /// ERROR: elements should not be enriched explicitly. 
        /// Instead the enrichment items should store which elements they interact with. 
        /// If the element needs to access the enrichment items it should do so through its nodes.
        /// </summary>
        public List<IEnrichmentItem2D> EnrichmentItems { get; }

        public XContinuumElement2D(IsoparametricElementType2D type, IReadOnlyList<XNode2D> nodes,
            IIntegrationStrategy2D<XContinuumElement2D> integrationStrategy)
        {
            type.CheckNodes(nodes);
            this.Nodes = nodes;
            this.elementType = type;
            this.EnrichmentItems = new List<IEnrichmentItem2D>();
            this.IntegrationStrategy = integrationStrategy;
        }

        public SymmetricMatrix2D BuildStandardStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D(StandardDofsCount);
            foreach (var gausspointMaterialPair in IntegrationStrategy.GetIntegrationPointsAndMaterials(this))
            {
                GaussPoint2D gaussPoint = gausspointMaterialPair.Key;
                IFiniteElementMaterial2D material = gausspointMaterialPair.Value;

                // Calculate the necessary quantities for the integration
                Matrix2D constitutive = material.CalculateConstitutiveMatrix();
                EvaluatedInterpolation2D evaluatedInterpolation =
                    elementType.Interpolation.EvaluateOnlyDerivativesAt(Nodes, gaussPoint);
                Matrix2D deformation = CalculateStandardDeformationMatrix(evaluatedInterpolation);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(material.Thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight); // Perhaps I shoul scale only the smallest matrix (constitutive) before the multiplications
                Debug.Assert(partial.Rows == StandardDofsCount);
                Debug.Assert(partial.Columns == StandardDofsCount);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(partial, stiffness);
            }
            return stiffness;
        }

        public void BuildEnrichedStiffnessMatrices(out Matrix2D stiffnessEnrichedStandard,
            out SymmetricMatrix2D stiffnessEnriched)
        {
            int standardDofsCount = StandardDofsCount;
            int artificialDofsCount = CountArtificialDofs();
            stiffnessEnrichedStandard = new Matrix2D(artificialDofsCount, standardDofsCount);
            stiffnessEnriched = new SymmetricMatrix2D(artificialDofsCount);

            foreach (var pair in IntegrationStrategy.GetIntegrationPointsAndMaterials(this))
            {
                GaussPoint2D gaussPoint = pair.Key;
                IFiniteElementMaterial2D material = pair.Value;

                // Calculate the necessary quantities for the integration
                Matrix2D constitutive = material.CalculateConstitutiveMatrix();
                EvaluatedInterpolation2D evaluatedInterpolation = Interpolation.EvaluateAt(Nodes, gaussPoint);
                Matrix2D Bstd = CalculateStandardDeformationMatrix(evaluatedInterpolation);
                Matrix2D Benr = CalculateEnrichedDeformationMatrix(artificialDofsCount,
                    gaussPoint, evaluatedInterpolation);

                // Contributions of this gauss point to the element stiffness matrices. 
                // Kee = SUM(Benr^T * E * Benr * dV), Kes = SUM(Benr^T * E * Bstd * dV)
                double dVolume = material.Thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight;
                Matrix2D transposeBenrTimesConstitutive = Benr.Transpose() * constitutive; // cache the result

                Matrix2D Kes = transposeBenrTimesConstitutive * Bstd;  // enriched-standard part
                Kes.Scale(dVolume); // TODO: Scale only the smallest matrix (constitutive) before the multiplications. Probably requires a copy of the constitutive matrix.
                MatrixUtilities.AddPartialToTotalMatrix(Kes, stiffnessEnrichedStandard);

                Matrix2D Kee = transposeBenrTimesConstitutive * Benr;  // enriched-enriched part
                Kee.Scale(dVolume);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(Kee, stiffnessEnriched);
            }
        }

        /// <summary>
        /// Calculates the deformation matrix B. Dimensions = 3x8.
        /// B is a linear transformation FROM the nodal values of the displacement field TO the the derivatives of
        /// the displacement field in respect to the cartesian axes (i.e. the stresses): {dU/dX} = [B] * {d} => 
        /// {u,x v,y u,y, v,x} = [... Bk ...] * {u1 v1 u2 v2 u3 v3 u4 v4}, where k = 1, ... nodesCount is a node and
        /// Bk = [dNk/dx 0; 0 dNk/dY; dNk/dy dNk/dx] (3x2)
        /// </summary>
        /// <param name="evaluatedInterpolation">The shape function derivatives calculated at a specific 
        ///     integration point</param>
        /// <returns></returns>
        public Matrix2D CalculateStandardDeformationMatrix(EvaluatedInterpolation2D evaluatedInterpolation)
        {
            var deformationMatrix = new Matrix2D(3, StandardDofsCount);
            for (int nodeIndex = 0; nodeIndex < Nodes.Count; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                Tuple<double, double> dNdX = evaluatedInterpolation.GetGlobalCartesianDerivativesOf(Nodes[nodeIndex]);

                deformationMatrix[0, col1] = dNdX.Item1;
                deformationMatrix[1, col2] = dNdX.Item2;
                deformationMatrix[2, col1] = dNdX.Item2;
                deformationMatrix[2, col2] = dNdX.Item1;
            }
            return deformationMatrix;
        }

        protected Matrix2D CalculateEnrichedDeformationMatrix(int artificialDofsCount,
            GaussPoint2D gaussPoint, EvaluatedInterpolation2D evaluatedInterpolation)
        {
            ICartesianPoint2D cartesianPoint = evaluatedInterpolation.TransformPointNaturalToGlobalCartesian(gaussPoint);
            var uniqueFunctions = new Dictionary<IEnrichmentFunction2D, EvaluatedFunction2D>();

            var deformationMatrix = new Matrix2D(3, artificialDofsCount);
            int currentColumn = 0;
            foreach (XNode2D node in Nodes)
            {
                double N = evaluatedInterpolation.GetValueOf(node);
                var dNdx = evaluatedInterpolation.GetGlobalCartesianDerivativesOf(node);

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

        /// <summary>
        /// The displacement field derivatives are a 2x2 matrix: gradientU[i,j] = dui/dj where i is the vector component 
        /// and j is the coordinate, w.r.t which the differentiation is done. The differentation coordinates and the
        /// vector components refer to the global cartesian system. 
        /// </summary>
        /// <param name="evaluatedInterpolation"></param>
        /// <param name="nodalDisplacementsX"></param>
        /// <param name="nodalDisplacementsY"></param>
        /// <returns></returns>
        protected DenseMatrix CalculateDisplacementFieldGradient(EvaluatedInterpolation2D evaluatedInterpolation, 
            double[] nodalDisplacements)
        {
            double[,] displacementGradient = new double[2, 2];
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double displacementX = nodalDisplacements[2 * nodeIdx];
                double displacementY = nodalDisplacements[2 * nodeIdx + 1];

                Tuple<double, double> shapeFunctionDerivatives = 
                    evaluatedInterpolation.GetGlobalCartesianDerivativesOf(Nodes[nodeIdx]);
                displacementGradient[0, 0] += shapeFunctionDerivatives.Item1 * displacementX;
                displacementGradient[0, 1] += shapeFunctionDerivatives.Item2 * displacementX;
                displacementGradient[1, 0] += shapeFunctionDerivatives.Item1 * displacementY;
                displacementGradient[1, 1] += shapeFunctionDerivatives.Item2 * displacementY;
            }
            return new DenseMatrix(displacementGradient);
        }

        // In a non linear problem I would also have to pass the new displacements or I would have to update the
        // material state elsewhere.
        protected Tensor2D CalculateStressTensor(DenseMatrix displacementFieldGradient, IFiniteElementMaterial2D material)
        {
            Matrix2D constitutive = material.CalculateConstitutiveMatrix(); 

            double strainXX = displacementFieldGradient[0, 0];
            double strainYY = displacementFieldGradient[1, 1];
            double strainXY = 0.5 * (displacementFieldGradient[0, 1] + displacementFieldGradient[1, 0]);

            // Should constitutive also be a tensor? Or  should I use matrices and vectors instead of tensors?
            double stressXX = constitutive[0, 0] * strainXX + constitutive[0, 1] * strainYY;
            double stressYY = constitutive[1, 0] * strainXX + constitutive[1, 1] * strainYY;
            double stressXY = constitutive[2, 2] * strainXY;

            return new Tensor2D(stressXX, stressYY, stressXY);
        }

        #region Dofs (perhaps all these should be delegated to element specific std and enr DofEnumerators)
        public int StandardDofsCount { get { return Nodes.Count * 2; } } // I could store it for efficency and update it when nodes change.

        public int CountArtificialDofs()
        {
            int count = 0;
            foreach (XNode2D node in Nodes) count += node.ArtificialDofsCount; // in all nodes or in enriched interpolation nodes?
            return count;
        }

        /// <summary>
        /// TODO: This should return readonly collections publicly.
        /// TODO: Perhaps this should be saved as a DOFEnumerator object (the dofs themselves would be created on  
        /// demand though). XElement will have a mutable one, while others will get a view. I could still use a  
        /// DOFEnumerator even if I do not save it. Transfering most of the code to the Enumerator class, also reduces  
        /// code duplication with the standard ContinuumElement2D
        /// </summary>
        /// <returns></returns>
        public IReadOnlyDictionary<Node2D, HashSet<StandardDOFType>> GetStandardNodalDOFTypes()
        {
            var nodalDOFTypes = new Dictionary<Node2D, HashSet<StandardDOFType>>(Nodes.Count);
            foreach (Node2D node in Nodes)
            {
                var dofTypesOfThisNode = new HashSet<StandardDOFType>();
                dofTypesOfThisNode.Add(StandardDOFType.X);
                dofTypesOfThisNode.Add(StandardDOFType.Y);
                nodalDOFTypes.Add(node, dofTypesOfThisNode);
            }
            return nodalDOFTypes;
        }
        #endregion
    }
}