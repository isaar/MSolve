using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Interpolation.ShapeFunctions;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    // Each element has its Interpolation object. Cons: Extra memory since there isn't any sharing (e.g. static 
    // interpolation classes). Pros: Each interpolation knows its nodes and perhaps node indices can be avoided.
    // The interpolation object creates a new InterpolationDerivatives at each integration point, which is promptly freed.
    class IsoparametricInterpolation2D
    {
        private readonly NaturalShapeFunctions2D shapeFunctions;

        public IReadOnlyList<Node2D> Nodes { get; } // Extra memory for a pointer only

        public IsoparametricInterpolation2D(IReadOnlyList<Node2D> nodes, NaturalShapeFunctions2D shapeFunctions)
        {
            if (nodes.Count != shapeFunctions.Count) // Perhaps needs more detailed checks.
            {
                throw new ArgumentException("The number of nodes and the shape functions do not match.");
            }
            this.Nodes = nodes;
            this.shapeFunctions = shapeFunctions;
        }

        public ShapeFunctionValues2D EvaluateAt(double xi, double eta)
        {
            return new ShapeFunctionValues2D(shapeFunctions.EvaluateAt(xi, eta), Nodes);
        }

        public ShapeFunctionDerivatives2D EvaluateDerivativesAt(double xi, double eta)
        {
            double[,] naturalDerivatives = shapeFunctions.EvaluateDerivativesAt(xi, eta);
            return new ShapeFunctionDerivatives2D(naturalDerivatives, new Jacobian2D(Nodes, naturalDerivatives));
        }
    }
}
