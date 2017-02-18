using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Elements
{
    class Quad4: IFiniteElement2DView
    {
        private const double DETERMINANT_TOLERANCE = 0.00000001;
        private const int DOFS_COUNT = 8;

        private readonly double halfLengthX;
        private readonly double halfLengthY;
        public IReadOnlyList<Node2D> Nodes { get; }
        private readonly IReadOnlyList<GaussPoint2D> gaussPoints;
        private readonly IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materials;

        //public IFiniteElementDOFEnumerator DOFEnumerator
        //{
        //    get { throw new NotImplementedException(); }
        //    set { throw new NotImplementedException(); }
        //}

        public Quad4(Node2D[] nodes, IFiniteElementMaterial2D material)
        {
            // TODO: Add checks here
            if ((nodes[0].X != nodes[3].X) || (nodes[1].X != nodes[2].X) 
                || (nodes[0].Y != nodes[1].Y) || (nodes[2].Y != nodes[3].Y))
            {
                throw new ArgumentException("The local cartsian system is not aligned to the global cartesian system or "
                    + "may not be even a rectangle. Use an isoparametric quad 4 instead or check the order of the nodes");
            }

            this.halfLengthX = 0.5 * (nodes[1].X - nodes[0].X);
            this.halfLengthY = 0.5 * (nodes[3].Y - nodes[0].Y);

            this.Nodes = new List<Node2D>(nodes);
            this.gaussPoints = FindGaussPoints();

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

        //public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        //{
        //    throw new NotImplementedException();
        //}

        //public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        //{
        //    throw new NotImplementedException();
        //}

        public SymmetricMatrix2D<double> BuildStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(DOFS_COUNT);
            foreach (var gaussPoint in gaussPoints) // TODO: remove the integration logic from the element class
            {
                // Calculate the necessary quantities for the integration
                Matrix2D<double> deformation = CalculateDeformationMatrix(gaussPoint);
                Matrix2D<double> constitutive = materials[gaussPoint].CalculateConstitutiveMatrix();
                double thickness = materials[gaussPoint].Thickness;

                // Gauss integration at this point
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation;
                partial.Scale(thickness * gaussPoint.Weight);
                Debug.Assert(partial.Rows == DOFS_COUNT);
                Debug.Assert(partial.Columns == DOFS_COUNT);
                AddPartialToSymmetricTotalMatrix(partial, stiffness);
            }
            return stiffness;
        }

        private static void AddPartialToSymmetricTotalMatrix(Matrix2D<double> partialMatrix,
            SymmetricMatrix2D<double> totalMatrix)
        {
            for (int row = 0; row < totalMatrix.Rows; ++row)
            {
                for (int col = row; col < totalMatrix.Columns; ++col)
                {
                    totalMatrix[row, col] += partialMatrix[row, col];
                }
            }
        }

        private IReadOnlyList<GaussPoint2D> FindGaussPoints()
        {
            var localGaussPoints = new List<GaussPoint2D>();
            foreach (var naturalGaussPoint in GaussQuadrature2D.Order2x2.GenerateIntegrationPoints())
            {
                double localX = naturalGaussPoint.X * halfLengthX;
                double localY = naturalGaussPoint.Y * halfLengthY;
                double localWeight = naturalGaussPoint.Weight * halfLengthX * halfLengthY;
                localGaussPoints.Add(new GaussPoint2D(localX, localY, localWeight));
            }
            return localGaussPoints;
        }

        /// <summary>
        /// Calculate the deformation matrix B (3x8).
        /// B is a linear transformation FROM the nodal values of the displacement field TO the strain vector: 
        /// {e} = [B] * {d} => {u,x v,y u,y+v,x} = [B] * {u1 v1 u2 v2 u3 v3 u4 v4}
        /// </summary>
        /// <param name="gaussPoint"></param>
        /// <returns>A 3x8 matrix</returns>
        private Matrix2D<double> CalculateDeformationMatrix(GaussPoint2D gaussPoint)
        {
            var B = new Matrix2D<double>(3, DOFS_COUNT);

            B[0, 0] = gaussPoint.Y - halfLengthY;
            B[0, 2] = -gaussPoint.Y + halfLengthY;
            B[0, 4] = gaussPoint.Y + halfLengthY;
            B[0, 6] = -gaussPoint.Y - halfLengthY;

            B[1, 1] = gaussPoint.X - halfLengthX;
            B[1, 3] = -gaussPoint.X - halfLengthX;
            B[1, 5] = gaussPoint.X + halfLengthX;
            B[1, 7] = -gaussPoint.X + halfLengthX;

            B[2, 0] = gaussPoint.X - halfLengthX;
            B[2, 1] = gaussPoint.Y - halfLengthY;
            B[2, 2] = -gaussPoint.X - halfLengthX;
            B[2, 3] = -gaussPoint.Y + halfLengthY;
            B[2, 4] = gaussPoint.X + halfLengthX;
            B[2, 5] = gaussPoint.Y + halfLengthY;
            B[2, 6] = -gaussPoint.X + halfLengthX;
            B[2, 7] = -gaussPoint.Y - halfLengthY;

            B.Scale(1.0 / (4 * halfLengthX * halfLengthY));
            return B;
        }
    }
}
