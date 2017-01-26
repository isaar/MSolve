using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Integration.ShapeFunctions
{
    /// <summary>
    /// The shape functions of a 4-node quadrilateral finite element.
    /// </summary>
    static class Quad4ShapeFunctions
    {
        private class BottomLeft : IShapeFunction2D
        {
            public double ValueAt(double xi, double eta)
            {
                return 0.25 * (1 - xi) * (1 - eta);
            }

            public double XiDerivativeAt(double xi, double eta)
            {
                return -0.25 * (1 - eta);
            }

            public double EtaDerivativeAt(double xi, double eta)
            {
                return -0.25 * (1 - xi);
            }
        }

        private class BottomRight : IShapeFunction2D
        {
            public double ValueAt(double xi, double eta)
            {
                return 0.25 * (1 + xi) * (1 - eta);
            }

            public double XiDerivativeAt(double xi, double eta)
            {
                return 0.25 * (1 - eta);
            }

            public double EtaDerivativeAt(double xi, double eta)
            {
                return -0.25 * (1 + xi);
            }
        }

        private class TopRight : IShapeFunction2D
        {
            public double ValueAt(double xi, double eta)
            {
                return 0.25 * (1 + xi) * (1 + eta);
            }

            public double XiDerivativeAt(double xi, double eta)
            {
                return 0.25 * (1 + eta);
            }

            public double EtaDerivativeAt(double xi, double eta)
            {
                return 0.25 * (1 + xi);
            }
        }

        private class TopLeft : IShapeFunction2D
        {
            public double ValueAt(double xi, double eta)
            {
                return 0.25 * (1 - xi) * (1 + eta);
            }

            public double XiDerivativeAt(double xi, double eta)
            {
                return -0.25 * (1 + eta);
            }

            public double EtaDerivativeAt(double xi, double eta)
            {
                return 0.25 * (1 - xi);
            }
        }

        /// <summary>
        /// Interpolates the values of (-1, -1)
        /// </summary>
        public static readonly IShapeFunction2D N1 = new BottomLeft();

        /// <summary>
        /// Interpolates the values of (1, -1)
        /// </summary>
        public static readonly IShapeFunction2D N2 = new BottomRight();

        /// <summary>
        /// Interpolates the values of (1, 1)
        /// </summary>
        public static readonly IShapeFunction2D N3 = new TopRight();

        /// <summary>
        /// Interpolates the values of (-1, 1)
        /// </summary>
        public static readonly IShapeFunction2D N4 = new TopLeft();

        public static IShapeFunction2D[] AllFunctions { get { return new IShapeFunction2D[] { N1, N2, N3, N4 }; } }

        public static double[] AllValuesAt(double xi, double eta)
        {
            double[] values = new double[4];
            IShapeFunction2D[] functions = AllFunctions;
            for (int i = 0; i < 4; ++i)
            {
                values[i] = functions[i].ValueAt(xi, eta);
            }
            return values;
        }

        /// <summary>
        /// Calculates the partial derivatives of all shape functions at the specified point. Each entry in the 
        /// collection is a (Ni_xi, Ni_eta) pair, where i is the index of the shape function, Ni_xi is the partial  
        /// derivative in respect to xi and Ni_eta is the partial derivative in respect to eta.
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        public static Tuple<double, double>[] AllDerivativesAt(double xi, double eta)
        {
            var derivatives = new Tuple<double, double>[4];
            IShapeFunction2D[] functions = AllFunctions;
            for (int i = 0; i < 4; ++i)
            {
                derivatives[i] = new Tuple<double, double>(
                    functions[i].XiDerivativeAt(xi, eta), functions[i].EtaDerivativeAt(xi, eta));
            }
            return derivatives;
        }
    }
}
