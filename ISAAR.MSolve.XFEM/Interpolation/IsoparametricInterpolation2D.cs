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

        public EvaluatedInterpolation2D EvaluateAt(IPoint2D naturalPoint)
        {
            double xi = naturalPoint.X;
            double eta = naturalPoint.Y;
            double[,] naturalDerivatives = shapeFunctions.EvaluateDerivativesAt(xi, eta);
            // TODO: perhaps check that the nodes match the values and derivatives
            return new EvaluatedInterpolation2D(Nodes, shapeFunctions.EvaluateAt(xi, eta), 
                naturalDerivatives, new Jacobian2D(Nodes, naturalDerivatives));
        }

        public EvaluatedInterpolation2D EvaluateOnlyDerivativesAt(IPoint2D naturalPoint)
        {
            double[,] naturalDerivatives = shapeFunctions.EvaluateDerivativesAt(naturalPoint.X, naturalPoint.Y);
            // TODO: perhaps check that the nodes match the values and derivatives
            return new EvaluatedInterpolation2D(Nodes, naturalDerivatives, new Jacobian2D(Nodes, naturalDerivatives));
        }
    }
}
