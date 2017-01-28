using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Integration.ShapeFunctions
{
    class Quad4ShapeFunctions
    {
        private readonly double halfLengthX;
        private readonly double halfLengthY;

        public Quad4ShapeFunctions(double halfLengthX, double halfLengthY)
        {
            this.halfLengthX = halfLengthX;
            this.halfLengthY = halfLengthY;
        }

        public double[] AllValuesAt(double x, double y)
        {
            double[] values = new double[4];
            values[0] = 0.25 * (1 - x / halfLengthX) * (1 - y / halfLengthY);
            values[1] = 0.25 * (1 + x / halfLengthX) * (1 - y / halfLengthY);
            values[2] = 0.25 * (1 + x / halfLengthX) * (1 + y / halfLengthY);
            values[3] = 0.25 * (1 - x / halfLengthX) * (1 + y / halfLengthY);
            return values;
        }

        /// <summary>
        /// Calculates the partial derivatives of all shape functions at the specified point. Each entry in the 
        /// collection is a (Ni_x, Ni_y) pair, where i is the index of the shape function, Ni_x is the partial  
        /// derivative in respect to x and Ni_y is the partial derivative in respect to y.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        public Tuple<double, double>[] AllDerivativesAt(double x, double y)
        {
            var derivatives = new Tuple<double, double>[4];
            derivatives[0] = new Tuple<double, double>(-0.25 / halfLengthX * (1 - y / halfLengthY), -0.25 / halfLengthY * (1 - x / halfLengthX));
            derivatives[1] = new Tuple<double, double>(0.25 / halfLengthX * (1 - y / halfLengthY), -0.25 / halfLengthY * (1 + x / halfLengthX));
            derivatives[2] = new Tuple<double, double>(0.25 / halfLengthX * (1 + y / halfLengthY), 0.25 / halfLengthY * (1 + x / halfLengthX));
            derivatives[3] = new Tuple<double, double>(-0.25 / halfLengthX * (1 + y / halfLengthY), 0.25 / halfLengthY * (1 - x / halfLengthX));
            return derivatives;
        }
    }
}
