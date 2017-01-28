using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Integration
{
    class Jacobian2D
    {
        private static readonly int dimension = 2;
        public double Determinant { get; }
        public Matrix2D<double> InverseJ { get; } // I need a Matrix view for this
        
        public Jacobian2D(IReadOnlyList<IPoint2D> nodes, Tuple<double, double>[] shapeFunctionNaturalDerivatives) 
        {
            //TODO: enforce nodes and gradients being the same length (e.g. using BiList)
            if (nodes.Count != shapeFunctionNaturalDerivatives.Length)
            {
                throw new ArgumentException("The nodal corrdinates do not match the shape function derivatives.");
            }

            // The original matrix is not stored. Only the inverse and the determinant
            Matrix2D<double> jacobianMatrix = CalculateJacobianMatrix(nodes, shapeFunctionNaturalDerivatives);
            Determinant = CalculateDeterminant(jacobianMatrix);
            InverseJ = CalculateInverseJacobian(jacobianMatrix);
        }

        /// <summary>
        /// Calculates the B1 matrix (the first part of the deformation matrix), which is a linear transformation FROM  
        /// the derivatives of the displacement field in respect to the natural axes TO the strain vector: 
        /// {e} = [B1] * {dU/dXi} => {u,x v,y u,y+v,x} = [B1] * {u,xi u,eta v,xi, v,eta}. The dimensions of B1 are 3x4.
        /// </summary>
        /// <returns>The B1 (3x4) matrix</returns>
        public Matrix2D<double> CalculateB1DeformationMatrix() //TODO: This probably only applies for linear problems
        {
            var B1 = new Matrix2D<double>(3, 4);
            B1[0, 0] = InverseJ[0, 0]; B1[0, 1] = InverseJ[0, 1];
            B1[1, 2] = InverseJ[1, 0]; B1[1, 3] = InverseJ[1, 1];
            B1[2, 0] = InverseJ[1, 0]; B1[2, 1] = InverseJ[1, 1];
            B1[2, 2] = InverseJ[0, 0]; B1[2, 3] = InverseJ[0, 1];
            return B1;
        }

        private static Matrix2D<double> CalculateJacobianMatrix(IReadOnlyList<IPoint2D> nodes, 
            Tuple<double, double>[] shapeFunctionNaturalDerivatives)
        {
            var J = new Matrix2D<double>(dimension, dimension);
            for (int n = 0; n < nodes.Count; ++n)
            {
                double x = nodes[n].X;
                double y = nodes[n].Y;
                double N_xi = shapeFunctionNaturalDerivatives[n].Item1;
                double N_eta = shapeFunctionNaturalDerivatives[n].Item2;

                J[0, 0] += N_xi * x;
                J[0, 1] += N_xi * y;
                J[1, 0] += N_eta * x;
                J[1, 1] += N_eta * y;
            }
            return J;
        }

        private static double CalculateDeterminant(Matrix2D<double> jacobianMatrix)
        {
            return jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[1, 0] * jacobianMatrix[0, 1];
        }

        private Matrix2D<double> CalculateInverseJacobian(Matrix2D<double> jacobianMatrix)
        {
            var invJ = new Matrix2D<double>(dimension, dimension);
            invJ[0, 0] = jacobianMatrix[1, 1] / Determinant;
            invJ[0, 1] = -jacobianMatrix[0, 1] / Determinant;
            invJ[1, 0] = -jacobianMatrix[1, 0] / Determinant;
            invJ[1, 1] = jacobianMatrix[0, 0] / Determinant;
            return invJ;
        }
    }
}
