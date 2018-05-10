using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Elements
{
    class ContinuumElement2D
    {
        private readonly IsoparametricElementType2D elementType;

        public IReadOnlyList<Node2D> Nodes { get; private set; }
        public IsoparametricInterpolation2D Interpolation { get { return elementType.Interpolation; } }
        public IStandardQuadrature2D StandardQuadrature { get { return elementType.StandardQuadrature; } }
        public IIntegrationStrategy2D<ContinuumElement2D> IntegrationStrategy { get; }
        public IMaterialField2D Material { get; }
        public int DofsCount { get { return Nodes.Count * 2; } } // I could store it for efficency and update it when nodes change.
        
        public ContinuumElement2D(IsoparametricElementType2D type, IReadOnlyList<Node2D> nodes,  
            IIntegrationStrategy2D<ContinuumElement2D> integrationStrategy, IMaterialField2D material)
        {
            type.CheckNodes(nodes);
            this.Nodes = nodes;
            this.elementType = type;
            this.IntegrationStrategy = integrationStrategy;
            this.Material = material;
        }

        public Matrix BuildStiffnessMatrix()
        {
            var stiffness = Matrix.CreateZero(DofsCount, DofsCount);
            foreach (GaussPoint2D gaussPoint in IntegrationStrategy.GenerateIntegrationPoints(this))
            {
                // Calculate the necessary quantities for the integration
                EvaluatedInterpolation2D evaluatedInterpolation =
                    elementType.Interpolation.EvaluateOnlyDerivativesAt(Nodes, gaussPoint);
                double thickness = Material.GetThicknessAt(gaussPoint, evaluatedInterpolation);
                Matrix constitutive = Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolation);
                Matrix deformation = CalculateDeformationMatrix(evaluatedInterpolation);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                double dVolume = thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight; // Perhaps I should scale only the smallest matrix (constitutive) before the multiplications 
                stiffness.AxpyIntoThis(partial, dVolume);
                Debug.Assert(partial.NumRows == DofsCount);
                Debug.Assert(partial.NumColumns == DofsCount);
                
            }
            return stiffness;
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
        public Matrix CalculateDeformationMatrix(EvaluatedInterpolation2D evaluatedInterpolation)
        {
            var deformationMatrix = Matrix.CreateZero(3, DofsCount);
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

        /// <summary>
        /// The displacement field derivatives are a 2x2 matrix: gradientU[i,j] = dui/dj where i is the vector component 
        /// and j is the coordinate, w.r.t which the differentiation is done. The differentation coordinates and the
        /// vector components refer to the global cartesian system. 
        /// </summary>
        /// <param name="evaluatedInterpolation"></param>
        /// <param name="nodalDisplacementsX"></param>
        /// <param name="nodalDisplacementsY"></param>
        /// <returns></returns>
        protected Matrix2by2 CalculateDisplacementFieldGradient(EvaluatedInterpolation2D evaluatedInterpolation,
            Vector nodalDisplacements)
        {
            var displacementGradient = Matrix2by2.CreateZero();
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                double displacementX = nodalDisplacements[2 * nodeIdx];
                double displacementY = nodalDisplacements[2 * nodeIdx + 1];

                Vector2 shapeFunctionDerivatives =
                    evaluatedInterpolation.GetGlobalCartesianDerivativesOf(Nodes[nodeIdx]);
                displacementGradient[0, 0] += shapeFunctionDerivatives[0] * displacementX;
                displacementGradient[0, 1] += shapeFunctionDerivatives[1] * displacementX;
                displacementGradient[1, 0] += shapeFunctionDerivatives[0] * displacementY;
                displacementGradient[1, 1] += shapeFunctionDerivatives[1] * displacementY;
            }
            return displacementGradient;
        }

        // In a non linear problem I would also have to pass the new displacements or I would have to update the
        // material state elsewhere.
        protected Tensor2D CalculateStressTensor(Matrix displacementFieldGradient, Matrix constitutive)
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

        // TODO: This should return readonly collections publicly, but the XElement class should still have access to 
        // the mutable collections
        // TODO: Perhaps this should be saved as a DofOrderer object. XElement will get a mutable one, while 
        // others will get a view. I could still use a DofOrderer even if I do not save it.
        public IReadOnlyDictionary<Node2D, HashSet<DisplacementDof>> GetNodalDofTypes()
        {
            var nodalDofTypes = new Dictionary<Node2D, HashSet<DisplacementDof>>(Nodes.Count);
            foreach (Node2D node in Nodes)
            {
                var dofTypesOfThisNode = new HashSet<DisplacementDof>();
                dofTypesOfThisNode.Add(DisplacementDof.X);
                dofTypesOfThisNode.Add(DisplacementDof.Y);
                nodalDofTypes.Add(node, dofTypesOfThisNode);
            }
            return nodalDofTypes;
        }
    }
}
