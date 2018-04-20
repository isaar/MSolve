using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    // TODO: The contnets of this class will almost always be accessed in for loops with the node index as a loop
    // counter, because there will generally be a need for the index's position in a matrix or there would be a list of
    // nodal values to interpolate. Using Node2D as key of the dictionaries with and requiring it from the client, 
    // introduces little extra safety, since the user will have to map the node index to Node2D. Performance wise it is
    // also a bad idea, since it requires to convert int -> Node2D (client) and then Dictionary lookup of Node2D (this 
    // class), instead of an array read. I guess it would be acceptable to use node indices to communicate with the 
    // client, IF that client is guaranteed to be the element. 
    class EvaluatedInterpolation2D
    {
        // TODO: these 2 dictionaries may be able to be optimized into 1. E.g. only 1 data structure 
        // or arrays with a node to index dictionary
        private readonly Dictionary<Node2D, double> nodesToValues;
        private readonly Dictionary<Node2D, Vector2> nodesToCartesianDerivatives;

        public Jacobian2D Jacobian { get; }

        public EvaluatedInterpolation2D(IReadOnlyList<Node2D> nodes, double[,] naturalDerivatives, Jacobian2D jacobian)
        {
            /// Any attempt at retrieving the not evaluated shape function values (through <see cref="GetValueOf"/> 
            /// will throw a NullReferenceException, which should be sufficient.
            nodesToValues = null;
             
            nodesToCartesianDerivatives = new Dictionary<Node2D, Vector2>(nodes.Count);
            for (int i = 0; i < nodes.Count; ++i)
            {
                nodesToCartesianDerivatives[nodes[i]] = jacobian.TransformNaturalDerivativesToCartesian(
                    Vector2.Create(naturalDerivatives[i, 0], naturalDerivatives[i, 1]));
            }
            this.Jacobian = jacobian;
        }

        public EvaluatedInterpolation2D(IReadOnlyList<Node2D> nodes, double[] shapeFunctionValues, 
            double[,] naturalDerivatives, Jacobian2D jacobian)
        {
            // TODO: Optimize the dictionaries. Could I provide comparers or sth that speeds up hashing?
            nodesToValues = new Dictionary<Node2D, double>(nodes.Count);
            nodesToCartesianDerivatives = new Dictionary<Node2D, Vector2>(nodes.Count);
            for (int i = 0; i < nodes.Count; ++i)
            {
                Node2D node = nodes[i];
                nodesToValues[node] = shapeFunctionValues[i];
                nodesToCartesianDerivatives[node] = jacobian.TransformNaturalDerivativesToCartesian(
                    Vector2.Create(naturalDerivatives[i, 0], naturalDerivatives[i, 1]));
            }
            this.Jacobian = jacobian;
        }

        public double GetValueOf(Node2D node)
        {
            return nodesToValues[node];
        }

        // TODO: It would be easier for the client to access directly dNi/dx and dNi/dy instead of a 
        // Tuple<double, double>. However it would require multiple node lookups. Otherwise, a software wide convention   
        // to use gradients as IReadOnlyList (immutable row arrays) or a dedicated (immutable class) can be enforced. 
        // TODO: Should the returned vector be immutable?
        public Vector2 GetGlobalCartesianDerivativesOf(Node2D node)
        {
            return nodesToCartesianDerivatives[node];
        }

        public ICartesianPoint2D TransformPointNaturalToGlobalCartesian(INaturalPoint2D naturalCoordinates)
        {
            double x = 0, y = 0;
            foreach (var entry in nodesToValues)
            {
                Node2D node = entry.Key;
                double val = entry.Value;
                x += val * node.X;
                y += val * node.Y;
            }
            return new CartesianPoint2D(x, y);
        }

        // TODO: Add methods: interpolate scalar, interpolate vector/point (there already is one, just make it similar 
        // to the other interpolations), interpolate scalar gradient, interpolate vector gradient. All these will need
        // nodal values as input. These should be passed in as immutable 2D and 2D lists for starters. Nodal 
        // dictionaries would require a lot of work for the client, but provide some extra security.
    }
}
