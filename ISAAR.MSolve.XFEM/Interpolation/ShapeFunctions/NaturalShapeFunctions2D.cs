using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Interpolation.ShapeFunctions
{
    /// <summary>
    /// Uses the Enum pattern.
    /// </summary>
    abstract class NaturalShapeFunctions2D
    {
        public static NaturalShapeFunctions2D Quad4 = new Quad4ShapeFunctions();

        // Prevents non-nested derived classes. What happens if the class gets too big?
        // i) Use an interface and many similar (and interchangeable) enum classes (e.g. 1 for quads, 1 for triangles).
        // ii) Relax the requirement that derived classes must be nested. 
        // iii) Keep the nested requirement but use partial class
        private NaturalShapeFunctions2D() { } 

        public abstract int Count { get; }
        public abstract double[] EvaluateAt(double xi, double eta);
        public abstract double[,] EvaluateDerivativesAt(double xi, double eta);

        private class Quad4ShapeFunctions: NaturalShapeFunctions2D
        {
            public override int Count { get { return 4; } }

            public override double[] EvaluateAt(double xi, double eta)
            {
                double[] values = new double[4];
                values[0] = 0.25 * (1 - xi) * (1 - eta);
                values[1] = 0.25 * (1 + xi) * (1 - eta);
                values[2] = 0.25 * (1 + xi) * (1 + eta);
                values[3] = 0.25 * (1 - xi) * (1 + eta);
                return values;
            }

            public override double[,] EvaluateDerivativesAt(double xi, double eta)
            {
                double[,] derivatives = new double[4, 2];
                derivatives[0, 0] = -0.25 * (1 - eta);
                derivatives[0, 1] = -0.25 * (1 - xi);
                derivatives[1, 0] = 0.25 * (1 - eta);
                derivatives[1, 1] = -0.25 * (1 + xi);
                derivatives[2, 0] = 0.25 * (1 + eta);
                derivatives[2, 1] = 0.25 * (1 + xi);
                derivatives[3, 0] = -0.25 * (1 + eta);
                derivatives[3, 1] = 0.25 * (1 - xi);
                return derivatives;
            }
        }
    }
}
