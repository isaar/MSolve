using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    class EvaluatedInterpolation2D
    {
        // TODO: these 2 dictionaries may be able to be optimized into 1. E.g. only 1 data structure 
        // or arrays with a node to index dictionary
        private readonly Dictionary<Node2D, double> nodesToValues;
        private readonly Dictionary<Node2D, Tuple<double, double>> nodesToCartesianDerivatives;

        public Jacobian2D Jacobian { get; }

        public EvaluatedInterpolation2D(IReadOnlyList<Node2D> nodes, double[,] naturalDerivatives, Jacobian2D jacobian)
        {
            /// Any attempt at retrieving the not evaluated shape function values (through <see cref="GetValueOf"/> 
            /// will throw a NullReferenceException, which should be sufficient.
            nodesToValues = null;
             
            nodesToCartesianDerivatives = new Dictionary<Node2D, Tuple<double, double>>(nodes.Count);
            for (int i = 0; i < nodes.Count; ++i)
            {
                nodesToCartesianDerivatives[nodes[i]] = jacobian.TransformNaturalDerivativesToCartesian(
                    naturalDerivatives[i, 0], naturalDerivatives[i, 1]);
            }
            this.Jacobian = jacobian;
        }

        public EvaluatedInterpolation2D(IReadOnlyList<Node2D> nodes, double[] shapeFunctionValues, 
            double[,] naturalDerivatives, Jacobian2D jacobian)
        {
            // TODO: Optimize the dictionaries. Could I provide comparers or sth that speeds up hashing?
            nodesToValues = new Dictionary<Node2D, double>(nodes.Count);
            nodesToCartesianDerivatives = new Dictionary<Node2D, Tuple<double, double>>(nodes.Count);
            for (int i = 0; i < nodes.Count; ++i)
            {
                Node2D node = nodes[i];
                nodesToValues[node] = shapeFunctionValues[i];
                nodesToCartesianDerivatives[node] = jacobian.TransformNaturalDerivativesToCartesian(
                    naturalDerivatives[i, 0], naturalDerivatives[i, 1]);
            }
            this.Jacobian = jacobian;
        }

        public double GetValueOf(Node2D node)
        {
            return nodesToValues[node];
        }

        public Tuple<double, double> GetCartesianDerivativesOf(Node2D node)
        {
            return nodesToCartesianDerivatives[node];
        }

        public IPoint2D TransformNaturalToCartesian(IPoint2D naturalCoordinates)
        {
            double x = 0, y = 0;
            foreach (var entry in nodesToValues)
            {
                Node2D node = entry.Key;
                double val = entry.Value;
                x += val * node.X;
                y += val * node.Y;
            }
            return new Point2D(x, y);
        }
    }
}
