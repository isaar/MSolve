using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Tensors;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.Discretization.FreedomDegrees;

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
    class XContinuumElement2D: ICell, IFiniteElement2DView<XNode2D, IsoparametricInterpolation2D>
    {
        public IsoparametricElementType2D ElementType { get; }

        public IReadOnlyList<ICartesianPoint2D> Vertices { get { return Nodes; } }

        /// <summary>
        /// All nodes are enriched for now.
        /// </summary>
        public IReadOnlyList<XNode2D> Nodes { get; }

        /// <summary>
        /// Common interpolation for standard and enriched nodes.
        /// </summary>
        public IsoparametricInterpolation2D Interpolation { get { return ElementType.Interpolation; } }

        public IStandardQuadrature2D StandardQuadrature { get { return ElementType.StandardQuadrature; } }
        public IIntegrationStrategy2D<XContinuumElement2D> IntegrationStrategy { get; }
        public IIntegrationStrategy2D<XContinuumElement2D> JintegralStrategy { get; }

        public IMaterialField2D Material { get; }

        /// <summary>
        /// ERROR: elements should not be enriched explicitly. 
        /// Instead the enrichment items should store which elements they interact with. 
        /// If the element needs to access the enrichment items it should do so through its nodes.
        /// Ok, but how would the integration strategy access the enrichment item?
        /// </summary>
        public List<IEnrichmentItem2D> EnrichmentItems { get; }

        public XContinuumElement2D(IsoparametricElementType2D type, IReadOnlyList<XNode2D> nodes,
            IMaterialField2D material, IIntegrationStrategy2D<XContinuumElement2D> integrationStrategy):
            this(type, nodes, material, integrationStrategy, integrationStrategy)
        {
            type.CheckNodes(nodes);
            this.Nodes = nodes;
            this.ElementType = type;
            this.EnrichmentItems = new List<IEnrichmentItem2D>();
            this.IntegrationStrategy = integrationStrategy;
            this.Material = material;
        }

        public XContinuumElement2D(IsoparametricElementType2D type, IReadOnlyList<XNode2D> nodes,
            IMaterialField2D material, IIntegrationStrategy2D<XContinuumElement2D> integrationStrategy,
            IIntegrationStrategy2D<XContinuumElement2D> jIntegralStrategy)
        {
            type.CheckNodes(nodes);
            this.Nodes = nodes;
            this.ElementType = type;
            this.EnrichmentItems = new List<IEnrichmentItem2D>();
            this.Material = material;
            this.IntegrationStrategy = integrationStrategy;
            this.JintegralStrategy = jIntegralStrategy;
        }

        // TODO: return a symmetric matrix
        public Matrix BuildStandardStiffnessMatrix()
        {
            var stiffness = Matrix.CreateZero(StandardDofsCount, StandardDofsCount);
            foreach (GaussPoint2D gaussPoint in IntegrationStrategy.GenerateIntegrationPoints(this))
            {
                // Calculate the necessary quantities for the integration
                EvaluatedInterpolation2D evaluatedInterpolation =
                    ElementType.Interpolation.EvaluateAt(Nodes, gaussPoint);
                double thickness = Material.GetThicknessAt(gaussPoint, evaluatedInterpolation);
                Matrix constitutive = Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolation);
                Matrix deformation = CalculateStandardDeformationMatrix(evaluatedInterpolation);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix partial = (deformation.MultiplyRight(constitutive, true, false)) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                Debug.Assert(partial.NumRows == StandardDofsCount);
                Debug.Assert(partial.NumColumns == StandardDofsCount);
                double dV = thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight; // Perhaps I should scale only the smallest matrix (constitutive) before the multiplications
                stiffness.AxpyIntoThis(partial, dV);
            }
            return stiffness;
        }

        public (Matrix stiffnessEnrichedEnriched, Matrix stiffnessEnrichedStandard) BuildEnrichedStiffnessMatricesLower()
        {
            int numStandardDofs = StandardDofsCount;
            int numEnrichedDofs = CountEnrichedDofs();
            var stiffnessEnrichedStandard = Matrix.CreateZero(numEnrichedDofs, numStandardDofs);
            var stiffnessEnrichedEnriched = Matrix.CreateZero(numEnrichedDofs, numEnrichedDofs);

            foreach (GaussPoint2D gaussPoint in IntegrationStrategy.GenerateIntegrationPoints(this))
            {
                // Calculate the necessary quantities for the integration
                EvaluatedInterpolation2D evaluatedInterpolation =
                    ElementType.Interpolation.EvaluateAt(Nodes, gaussPoint);
                double thickness = Material.GetThicknessAt(gaussPoint, evaluatedInterpolation);
                Matrix constitutive = Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolation);
                Matrix Bstd = CalculateStandardDeformationMatrix(evaluatedInterpolation);
                Matrix Benr = CalculateEnrichedDeformationMatrix(numEnrichedDofs,
                    gaussPoint, evaluatedInterpolation);

                // Contributions of this gauss point to the element stiffness matrices. 
                // Kee = SUM(Benr^T * E * Benr * dV), Kes = SUM(Benr^T * E * Bstd * dV)
                // TODO: Scale only the smallest matrix (constitutive) before the multiplications. Probably requires a copy of the constitutive matrix.
                double dVolume = thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight;
                Matrix transposeBenrTimesConstitutive = Benr.MultiplyRight(constitutive, true, false); // cache the result

                Matrix Kes = transposeBenrTimesConstitutive * Bstd;  // enriched-standard part
                stiffnessEnrichedStandard.AxpyIntoThis(Kes, dVolume);

                Matrix Kee = transposeBenrTimesConstitutive * Benr;  // enriched-enriched part
                stiffnessEnrichedEnriched.AxpyIntoThis(Kee, dVolume);
            }

            return (stiffnessEnrichedEnriched, stiffnessEnrichedStandard);
        }

        public (Matrix stiffnessEnrichedEnriched, Matrix stiffnessStandardEnriched) BuildEnrichedStiffnessMatricesUpper()
        {
            int numStandardDofs = StandardDofsCount;
            int numEnrichedDofs = CountEnrichedDofs();
            var stiffnessStandardEnriched = Matrix.CreateZero(numStandardDofs, numEnrichedDofs);
            var stiffnessEnrichedEnriched = Matrix.CreateZero(numEnrichedDofs, numEnrichedDofs);

            foreach (GaussPoint2D gaussPoint in IntegrationStrategy.GenerateIntegrationPoints(this))
            {
                // Calculate the necessary quantities for the integration
                EvaluatedInterpolation2D evaluatedInterpolation =
                    ElementType.Interpolation.EvaluateAt(Nodes, gaussPoint);
                double thickness = Material.GetThicknessAt(gaussPoint, evaluatedInterpolation);
                Matrix constitutive = Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolation);
                Matrix Bstd = CalculateStandardDeformationMatrix(evaluatedInterpolation);
                Matrix Benr = CalculateEnrichedDeformationMatrix(numEnrichedDofs,
                    gaussPoint, evaluatedInterpolation);

                // Contributions of this gauss point to the element stiffness matrices. 
                // Kee = SUM(Benr^T * E * Benr * dV), Kse = SUM(Bstd^T * E * Benr * dV)
                // TODO: Scale only the smallest matrix (constitutive) before the multiplications. Probably requires a copy of the constitutive matrix.
                double dVolume = thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight;
                Matrix constitutiveTimesBenr = constitutive * Benr; // cache the result

                Matrix Kse = Bstd.MultiplyRight(constitutiveTimesBenr, true, false);  // enriched-standard part
                stiffnessStandardEnriched.AxpyIntoThis(Kse, dVolume);

                Matrix Kee = Benr.MultiplyRight(constitutiveTimesBenr, true, false);  // enriched-enriched part
                stiffnessEnrichedEnriched.AxpyIntoThis(Kee, dVolume);
            }

            return (stiffnessEnrichedEnriched, stiffnessStandardEnriched);
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
        public Matrix CalculateStandardDeformationMatrix(EvaluatedInterpolation2D evaluatedInterpolation)
        {
            var deformationMatrix = Matrix.CreateZero(3, StandardDofsCount);
            for (int nodeIndex = 0; nodeIndex < Nodes.Count; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                Vector2 dNdX = evaluatedInterpolation.GetGlobalCartesianDerivativesOf(Nodes[nodeIndex]);

                deformationMatrix[0, col1] = dNdX[0];
                deformationMatrix[1, col2] = dNdX[1];
                deformationMatrix[2, col1] = dNdX[1];
                deformationMatrix[2, col2] = dNdX[0];
            }
            return deformationMatrix;
        }

        // TODO: the argument asrtificialDofsCount was added when this method was private and only called by 
        // BuildStiffnessMatrix() that already counted the dofs. Since it is now used by other modules 
        // (J-integral, output), it would be better to obscure it, at the cost of recounting the dofs in some cases.
        public Matrix CalculateEnrichedDeformationMatrix(int artificialDofsCount,
            INaturalPoint2D gaussPoint, EvaluatedInterpolation2D evaluatedInterpolation)
        {
            //ICartesianPoint2D cartesianPoint = evaluatedInterpolation.TransformPointNaturalToGlobalCartesian(gaussPoint);
            var uniqueEnrichments = new Dictionary<IEnrichmentItem2D, EvaluatedFunction2D[]>();

            var deformationMatrix = Matrix.CreateZero(3, artificialDofsCount);
            int currentColumn = 0;
            foreach (XNode2D node in Nodes)
            {
                double N = evaluatedInterpolation.GetValueOf(node);
                var dNdx = evaluatedInterpolation.GetGlobalCartesianDerivativesOf(node);

                foreach (var enrichment in node.EnrichmentItems)
                {
                    IEnrichmentItem2D enrichmentItem = enrichment.Key;
                    double[] nodalEnrichmentValues = enrichment.Value;

                    // The enrichment function probably has been evaluated when processing a previous node. Avoid reevaluation.
                    EvaluatedFunction2D[] evaluatedEnrichments;
                    if (!(uniqueEnrichments.TryGetValue(enrichmentItem, out evaluatedEnrichments)))
                    {
                        evaluatedEnrichments = enrichmentItem.EvaluateAllAt(gaussPoint, this, evaluatedInterpolation);
                        uniqueEnrichments[enrichmentItem] = evaluatedEnrichments;
                    }

                    for (int i = 0; i < evaluatedEnrichments.Length; ++i)
                    {
                        // For each node and with all derivatives w.r.t. cartesian coordinates, the enrichment derivatives 
                        // are: Bx = enrN,x = N,x(x,y) * [H(x,y) - H(node)] + N(x,y) * H,x(x,y), where H is the enrichment 
                        // function
                        double Bx = dNdx[0] * (evaluatedEnrichments[i].Value - nodalEnrichmentValues[i])
                            + N * evaluatedEnrichments[i].CartesianDerivatives[0];
                        double By = dNdx[1] * (evaluatedEnrichments[i].Value - nodalEnrichmentValues[i])
                            + N * evaluatedEnrichments[i].CartesianDerivatives[1];

                        // This depends on the convention: node major or enrichment major. The following is node major.
                        int col1 = currentColumn++;
                        int col2 = currentColumn++;

                        deformationMatrix[0, col1] = Bx;
                        deformationMatrix[1, col2] = By;
                        deformationMatrix[2, col1] = By;
                        deformationMatrix[2, col2] = Bx;
                    }
                }
            }
            Debug.Assert(currentColumn == artificialDofsCount);
            return deformationMatrix;
        }

        /// <summary>
        /// This only works for points that do not lie on the crack interface. As such it is safe to pass GPs only
        /// </summary>
        /// <param name="gaussPoint"></param>
        /// <param name="evaluatedInterpolation"></param>
        /// <param name="standardNodalDisplacements"></param>
        /// <param name="enrichedNodalDisplacements"></param>
        /// <returns></returns>
        public Vector2 CalculateDisplacementField(INaturalPoint2D gaussPoint, EvaluatedInterpolation2D evaluatedInterpolation,
            Vector standardNodalDisplacements, Vector enrichedNodalDisplacements)
        {
            #region debug
            double tol = 1e-6;
            if ((Math.Abs(gaussPoint.Xi) <= tol) || (Math.Abs(gaussPoint.Eta) <= tol))
            {
                Console.WriteLine("Found an intersection point that isn't a node.");
            }
            #endregion
            var displacements = new double[2];

            // Standard contributions
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double shapeFunction = evaluatedInterpolation.GetValueOf(Nodes[nodeIdx]);
                displacements[0] += shapeFunction * standardNodalDisplacements[2 * nodeIdx];
                displacements[1] += shapeFunction * standardNodalDisplacements[2 * nodeIdx + 1];
            }

            // Enriched contributions
            //TODO: this should be taken as input, so that it is only computed once if it is needed for displacements, strains, 
            //      stresses, just like the evaluated interpolation.
            IReadOnlyDictionary<IEnrichmentItem2D, EvaluatedFunction2D[]> evalEnrichments =
                EvaluateEnrichments(gaussPoint, evaluatedInterpolation); 
            int dof = 0;
            foreach (XNode2D node in Nodes)
            {
                double shapeFunction = evaluatedInterpolation.GetValueOf(node);
                foreach (var nodalEnrichment in node.EnrichmentItems)
                {
                    EvaluatedFunction2D[] currentEvalEnrichments = evalEnrichments[nodalEnrichment.Key];
                    for (int e = 0; e < currentEvalEnrichments.Length; ++e)
                    {
                        double basisFunction = shapeFunction * (currentEvalEnrichments[e].Value - nodalEnrichment.Value[e]);
                        //double basisFunction = shapeFunction * currentEvalEnrichments[e].Value; // for debugging
                        displacements[0] += basisFunction * enrichedNodalDisplacements[dof++];
                        displacements[1] += basisFunction * enrichedNodalDisplacements[dof++];
                    }
                }
            }

            return Vector2.CreateFromArray(displacements);
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
        public Matrix2by2 CalculateDisplacementFieldGradient(INaturalPoint2D gaussPoint, 
            EvaluatedInterpolation2D evaluatedInterpolation, Vector standardNodalDisplacements,
            Vector enrichedNodalDisplacements) //TODO: this must only allow evaluations at Gauss points. It doesn't work for points on the crack interface
        {
            var displacementGradient = Matrix2by2.CreateZero();

            // Standard contributions
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double displacementX = standardNodalDisplacements[2 * nodeIdx];
                double displacementY = standardNodalDisplacements[2 * nodeIdx + 1];

                Vector2 shapeFunctionDerivatives = 
                    evaluatedInterpolation.GetGlobalCartesianDerivativesOf(Nodes[nodeIdx]);
                displacementGradient[0, 0] += shapeFunctionDerivatives[0] * displacementX;
                displacementGradient[0, 1] += shapeFunctionDerivatives[1] * displacementX;
                displacementGradient[1, 0] += shapeFunctionDerivatives[0] * displacementY;
                displacementGradient[1, 1] += shapeFunctionDerivatives[1] * displacementY;
            }

            // Enriched contributions. TODO: Extract the common steps with building B into a separate method 
            IReadOnlyDictionary<IEnrichmentItem2D, EvaluatedFunction2D[]> evalEnrichments = 
                EvaluateEnrichments(gaussPoint, evaluatedInterpolation);
            int dof = 0;
            foreach (XNode2D node in Nodes)
            {
                double N = evaluatedInterpolation.GetValueOf(node);
                Vector2 gradN = evaluatedInterpolation.GetGlobalCartesianDerivativesOf(node);

                foreach (var nodalEnrichment in node.EnrichmentItems)
                {
                    EvaluatedFunction2D[] currentEvalEnrichments = evalEnrichments[nodalEnrichment.Key];
                    for (int e = 0; e < currentEvalEnrichments.Length; ++e)
                    {
                        double psi = currentEvalEnrichments[e].Value;
                        Vector2 gradPsi = currentEvalEnrichments[e].CartesianDerivatives;
                        double deltaPsi = psi - nodalEnrichment.Value[e];

                        double Bx = gradN[0] * deltaPsi + N * gradPsi[0];
                        double By = gradN[1] * deltaPsi + N * gradPsi[1];

                        double enrDisplacementX = enrichedNodalDisplacements[dof++];
                        double enrDisplacementY = enrichedNodalDisplacements[dof++];

                        displacementGradient[0, 0] += Bx * enrDisplacementX;
                        displacementGradient[0, 1] += By * enrDisplacementX;
                        displacementGradient[1, 0] += Bx * enrDisplacementY;
                        displacementGradient[1, 1] += By * enrDisplacementY;
                    }
                }
            }

            return displacementGradient;
        }

        // In a non linear problem I would also have to pass the new displacements or I would have to update the
        // material state elsewhere.
        public Tensor2D CalculateStressTensor(Matrix2by2 displacementFieldGradient, Matrix constitutive)
        {
            double strainXX = displacementFieldGradient[0, 0];
            double strainYY = displacementFieldGradient[1, 1];
            double strainXYtimes2 = displacementFieldGradient[0, 1] + displacementFieldGradient[1, 0];

            // Should constitutive also be a tensor? Or  should I use matrices and vectors instead of tensors?
            double stressXX = constitutive[0, 0] * strainXX + constitutive[0, 1] * strainYY;
            double stressYY = constitutive[1, 0] * strainXX + constitutive[1, 1] * strainYY;
            double stressXY = constitutive[2, 2] * strainXYtimes2;

            return new Tensor2D(stressXX, stressYY, stressXY);
        }

        #region Dofs (perhaps all these should be delegated to element specific std and enr DofOrderers)
        public int CountEnrichedDofs()
        {
            int count = 0;
            foreach (XNode2D node in Nodes) count += node.EnrichedDofsCount; // in all nodes or in enriched interpolation nodes?
            return count;
        }

        public int StandardDofsCount { get { return Nodes.Count * 2; } } // I could store it for efficency and update it when nodes change.

        public FreedomDegrees.Ordering.DofTable<EnrichedDof> GetEnrichedDofs()
        {
            var elementDofs = new FreedomDegrees.Ordering.DofTable<EnrichedDof>();
            int dofCounter = 0;
            foreach (XNode2D node in Nodes)
            {
                foreach (var enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (var enrichedDof in enrichment.Dofs) // there are different dofs for x and y axes
                    {
                        elementDofs[node, enrichedDof] = dofCounter++;
                    }
                }
            }
            return elementDofs;
        }

        /// <summary>
        /// TODO: Perhaps this should be saved as a DofOrderer object (the dofs themselves would be created on  
        /// demand though). XElement will have a mutable one, while others will get a view. I could still use a  
        /// DofOrderer even if I do not save it. Transfering most of the code to the Enumerator class, also reduces  
        /// code duplication with the standard ContinuumElement2D
        /// </summary>
        /// <returns></returns>
        public FreedomDegrees.Ordering.DofTable<DisplacementDof> GetStandardDofs()
        {
            var elementDofs = new FreedomDegrees.Ordering.DofTable<DisplacementDof>();
            int dofCounter = 0;
            foreach (XNode2D node in Nodes)
            {
                elementDofs[node, DisplacementDof.X] = dofCounter++;
                elementDofs[node, DisplacementDof.Y] = dofCounter++;
            }
            return elementDofs;
        }
        #endregion
        
        private IReadOnlyDictionary<IEnrichmentItem2D, EvaluatedFunction2D[]> EvaluateEnrichments(
            INaturalPoint2D gaussPoint, EvaluatedInterpolation2D evaluatedInterpolation)
        {
            var cachedEvalEnrichments = new Dictionary<IEnrichmentItem2D, EvaluatedFunction2D[]>();
            foreach (XNode2D node in Nodes)
            {
                foreach (var enrichment in node.EnrichmentItems)
                {
                    IEnrichmentItem2D enrichmentItem = enrichment.Key;
                    double[] nodalEnrichmentValues = enrichment.Value;

                    // The enrichment function probably has been evaluated when processing a previous node. Avoid reevaluation.
                    if (!(cachedEvalEnrichments.TryGetValue(enrichmentItem, out EvaluatedFunction2D[] evaluatedEnrichments)))
                    {
                        evaluatedEnrichments = enrichmentItem.EvaluateAllAt(gaussPoint, this, evaluatedInterpolation);
                        cachedEvalEnrichments[enrichmentItem] = evaluatedEnrichments;
                    }
                }
            }
            return cachedEvalEnrichments;
        }
    }
}