using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
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
        protected static double determinantTolerance = 0.00000001;

        private readonly IReadOnlyList<Node2D> nodes;
        private readonly IReadOnlyList<GaussPoint2D> gaussPoints;
        private readonly IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materials;

        //public IFiniteElementDOFEnumerator DOFEnumerator
        //{
        //    get { throw new NotImplementedException(); }
        //    set { throw new NotImplementedException(); }
        //}

        public IsoparametricQuad4(Node2D[] nodes, IFiniteElementMaterial2D material)
        {
            // TODO: Add checks here: order of nodes
            this.nodes = new List<Node2D>(nodes);
            this.gaussPoints = IntegrationRule2D.Order2x2.Points;

            var materialsDict = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in gaussPoints)
            {
                materialsDict[point] = material.Clone();
            }
            this.materials = materialsDict;
        }

        //public IList<IList<DOFType>> GetElementDOFTypes(Element element)
        //{
        //    throw new NotImplementedException();
        //}

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public SymmetricMatrix2D<double> BuildStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(8);
            foreach (var gaussPoint in gaussPoints) // TODO: remove the integration logic from the element class
            {
                // Calculate the necessary quantities for the integration
                Tuple<double, double>[] shapeFunctionDerivatives = 
                    Quad4ShapeFunctions.AllDerivativesAt(gaussPoint.X, gaussPoint.Y);
                Jacobian2D jacobian = new Jacobian2D(nodes, shapeFunctionDerivatives);
                if (jacobian.Determinant < determinantTolerance)
                {
                    throw new ArgumentException(String.Format(
                        "Jacobian determinant is negative or under tolerance ({0} < {1}). Check the order of nodes or the element geometry.",
                        jacobian.Determinant, determinantTolerance));
                }

                Matrix2D<double> deformation = CalculateDeformationMatrix(shapeFunctionDerivatives, jacobian);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Gauss integration at this point
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(thickness * jacobian.Determinant * gaussPoint.Weight);
                Debug.Assert(partial.Rows == 8);
                Debug.Assert(partial.Columns == 8);
                AddPartialMatrixToStiffness(stiffness, partial);
            }
            return stiffness;
        }

        private static void AddPartialMatrixToStiffness(SymmetricMatrix2D<double> stiffnessMatrix, 
            Matrix2D<double> partialMatrix)
        {
            for (int row = 0; row < stiffnessMatrix.Rows; ++row)
            {
                for (int col = row; col < stiffnessMatrix.Columns; ++col)
                {
                    stiffnessMatrix[row, col] += partialMatrix[row, col];
                }
            }
        }

        private Matrix2D<double> CalculateDeformationMatrix(Tuple<double, double>[] shapeFunctionDerivatives, 
            Jacobian2D jacobian)
        {
            //Calculate B2 deformation matrix. Dimensions = 4x8.
            //B2 is a linear transformation FROM the nodal values of the displacement field TO the the derivatives of
            //the displacement field in respect to the natural axes: {dU/dXi} = [B2] * {d} => 
            //{u,xi u,eta v,xi, v,eta} = [B2] * {u1 v1 u2 v2 u3 v3 u4 v4}
            var B2 = new Matrix2D<double>(4, 8);
            for (int node = 0; node < 4; ++node)
            {
                int col1 = 2 * node;
                int col2 = 2 * node + 1;
                double Nxi = shapeFunctionDerivatives[node].Item1;
                double Neta = shapeFunctionDerivatives[node].Item2;

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
