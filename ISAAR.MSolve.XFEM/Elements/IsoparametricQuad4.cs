using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Integration.GaussPoints;
using ISAAR.MSolve.XFEM.Integration.ShapeFunctions;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    class IsoparametricQuad4: IFiniteElement2D
    {
        private const int DOFS_COUNT = 8;

        private readonly IReadOnlyList<Node2D> nodes;
        private readonly IReadOnlyList<GaussPoint2D> gaussPoints;
        private readonly IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materials;

        public IsoparametricQuad4(Node2D[] nodes, IFiniteElementMaterial2D material)
        {
            // TODO: Add checks here: order of nodes
            this.nodes = new List<Node2D>(nodes);
            this.gaussPoints = IntegrationRule2D.Order2x2.Points; // TODO: remove the integration logic from the element class

            var materialsDict = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in gaussPoints)
            {
                materialsDict[point] = material.Clone();
            }
            this.materials = materialsDict;
        }

        public SymmetricMatrix2D<double> BuildStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(DOFS_COUNT);
            foreach (var gaussPoint in gaussPoints)
            {
                // Calculate the necessary quantities for the integration
                ShapeFunctionDerivatives2D shapeFunctionDerivatives = 
                    IsoparametricQuad4ShapeFunctions.AllDerivativesAt(gaussPoint.X, gaussPoint.Y);
                Jacobian2D jacobian = new Jacobian2D(nodes, shapeFunctionDerivatives);

                Matrix2D<double> deformation = CalculateDeformationMatrix(shapeFunctionDerivatives, jacobian);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Gauss integration at this point
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(thickness * jacobian.Determinant * gaussPoint.Weight);
                Debug.Assert(partial.Rows == DOFS_COUNT);
                Debug.Assert(partial.Columns == DOFS_COUNT);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(partial, stiffness);
            }
            return stiffness;
        }

        private Matrix2D<double> CalculateDeformationMatrix(ShapeFunctionDerivatives2D shapeFunctionDerivatives, 
            Jacobian2D jacobian)
        {
            //Calculate B2 deformation matrix. Dimensions = 4x8.
            //B2 is a linear transformation FROM the nodal values of the displacement field TO the the derivatives of
            //the displacement field in respect to the natural axes: {dU/dXi} = [B2] * {d} => 
            //{u,xi u,eta v,xi, v,eta} = [B2] * {u1 v1 u2 v2 u3 v3 u4 v4}
            var B2 = new Matrix2D<double>(4, DOFS_COUNT);
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
    }
}
