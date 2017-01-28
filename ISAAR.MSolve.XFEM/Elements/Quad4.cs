using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Integration.GaussPoints;
using ISAAR.MSolve.XFEM.Integration.ShapeFunctions;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Elements
{
    class Quad4: IFiniteElement2D
    {
        private const double DETERMINANT_TOLERANCE = 0.00000001;
        private const int DOFS_COUNT = 8;

        private readonly double halfLengthX;
        private readonly double halfLengthY;
        private readonly double centroidX;
        private readonly double centroidY;
        private readonly IReadOnlyList<Node2D> nodes;
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
            this.centroidX = 0.5 * (nodes[0].X + nodes[1].X);
            this.halfLengthY = 0.5 * (nodes[3].Y - nodes[0].Y);
            this.centroidY = 0.5 * (nodes[0].Y + nodes[3].Y);

            this.nodes = new List<Node2D>(nodes);
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
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
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
            foreach (var naturalGaussPoint in IntegrationRule2D.Order2x2.Points)
            {
                double localX = naturalGaussPoint.X * halfLengthX + centroidX;
                double localY = naturalGaussPoint.Y * halfLengthY + centroidY;
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
            var shapeFunctions = new Quad4ShapeFunctions(halfLengthX, halfLengthY);
            Tuple<double, double>[] shapeFunctionDerivatives = 
                shapeFunctions.AllDerivativesAt(gaussPoint.X, gaussPoint.Y);

            var B = new Matrix2D<double>(3, DOFS_COUNT);
            for (int node = 0; node < 4; ++node)
            {
                int col1 = 2 * node;
                int col2 = 2 * node + 1;
                double Nx = shapeFunctionDerivatives[node].Item1;
                double Ny = shapeFunctionDerivatives[node].Item2;

                B[0, col1] = Nx;
                B[1, col2] = Ny;
                B[2, col1] = Ny;
                B[2, col2] = Nx;
            }
            return B;
        }
    }
}
