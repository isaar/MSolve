using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    abstract class ContinuumElement2D
    {
        public IReadOnlyList<Node2D> Nodes { get; private set; }
        public IIntegrationStrategy2D IntegrationStrategy { get; }

        public int DofsCount { get { return Nodes.Count * 2; } } // I could store it for efficency and update it when nodes change.
        
        public IsoparametricInterpolation2D Interpolation { get; }

        /// <summary>
        /// TODO: Create a standard integration rule interface that guarantees gauss points that are
        /// i) immutable, ii) precached for fast generation, iii) stored globally for all elements
        /// TODO: Should this be stored as a field? It is a static property of the concrete (aka derived) class...
        /// </summary>
        public IIntegrationRule2D MinimumIntegrationRule { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="nodes">The caller is responsibile for checking their suitability 
        ///     (at least their number).</param>
        /// <param name="integrationStrategyFactory"></param>
        protected ContinuumElement2D(IReadOnlyList<Node2D> nodes, IsoparametricInterpolation2D interpolation, 
            IIntegrationRule2D minimumIntegrationRule, IIntegrationStrategyFactory2D integrationStrategyFactory)
        {
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.MinimumIntegrationRule = minimumIntegrationRule;

            /// WARNING: this will probably try to access the members of <see cref="ContinuumElement2D"/>. 
            /// Thus they must be ready. However it is easy to have members that are initialized in the constructor of 
            /// the derived class, but are still uninitialized when they are accessed
            this.IntegrationStrategy = integrationStrategyFactory.CreateStrategy(this); 
        }

        public SymmetricMatrix2D<double> BuildStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(DofsCount);
            foreach (var gausspointMaterialPair in IntegrationStrategy.GetIntegrationPointsAndMaterials())
            {
                GaussPoint2D gaussPoint = gausspointMaterialPair.Item1;
                IFiniteElementMaterial2D material = gausspointMaterialPair.Item2;

                // Calculate the necessary quantities for the integration
                Matrix2D<double> constitutive = material.CalculateConstitutiveMatrix();
                EvaluatedInterpolation2D evaluatedInterpolation =
                    Interpolation.EvaluateOnlyDerivativesAt(Nodes, gaussPoint);
                Matrix2D<double> deformation = CalculateDeformationMatrix(evaluatedInterpolation);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(material.Thickness * evaluatedInterpolation.Jacobian.Determinant * gaussPoint.Weight); // Perhaps I shoul scale only the smallest matrix (constitutive) before the multiplications
                Debug.Assert(partial.Rows == DofsCount);
                Debug.Assert(partial.Columns == DofsCount);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(partial, stiffness);
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
        public Matrix2D<double> CalculateDeformationMatrix(EvaluatedInterpolation2D evaluatedInterpolation)
        {
            var deformationMatrix = new Matrix2D<double>(3, DofsCount);
            for (int nodeIndex = 0; nodeIndex < Nodes.Count; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                Tuple<double, double> dNdX = evaluatedInterpolation.GetCartesianDerivativesOf(Nodes[nodeIndex]);

                deformationMatrix[0, col1] = dNdX.Item1;
                deformationMatrix[1, col2] = dNdX.Item2;
                deformationMatrix[2, col1] = dNdX.Item2;
                deformationMatrix[2, col2] = dNdX.Item1;
            }
            return deformationMatrix;
        }

        // TODO: This should return readonly collections publicly, but the XElement class should still have access to 
        // the mutable collections
        // TODO: Perhaps this should be saved as a DOFEnumerator object. XElement will get a mutable one, while 
        // others will get a view. I could still use a DOFEnumerator even if I do not save it.
        public IReadOnlyDictionary<Node2D, HashSet<StandardDOFType>> GetNodalDOFTypes()
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
    }
}
