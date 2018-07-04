using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Interpolation.GaussPointSystems
{
    // TODO: Use the existing interoplation classes. 
    // Would be ideal to use the existing quadrature too, but the order is different.
    class Quad4GPSystem: IGaussPointSystem
    {
        private const int numGPs = 4;
        private const int numNodes = 4;

        // The order is important
        private static readonly GaussPoint2D[] gaussPoints = new GaussPoint2D[]
        {
            new GaussPoint2D(-1/Math.Sqrt(3), -1/Math.Sqrt(3), 1.0),
            new GaussPoint2D(1/Math.Sqrt(3), -1/Math.Sqrt(3), 1.0),
            new GaussPoint2D(-1/Math.Sqrt(3), 1/Math.Sqrt(3), 1.0),
            new GaussPoint2D(1/Math.Sqrt(3), 1/Math.Sqrt(3), 1.0)
        };

        private static readonly Func<double, double, double>[] shapeFunctions = new Func<double, double, double>[]
        {
            (r, s) => { return 0.25 * (1 - r) * (1 - s); },
            (r, s) => { return 0.25 * (1 + r) * (1 - s); },
            (r, s) => { return 0.25 * (1 - r) * (1 + s); },
            (r, s) => { return 0.25 * (1 + r) * (1 + s); }            
        };

        private static readonly double[,] shapeFunctionsAtNodes;

        static Quad4GPSystem()
        {
            // Coordinates of nodes in the GP system: r=sqrt(3)*xi, s=sqrt(3)*eta
            double[,] nodeCoordinates = new double[,]
            {
                { -Math.Sqrt(3), -Math.Sqrt(3) },
                { Math.Sqrt(3), -Math.Sqrt(3) },
                { Math.Sqrt(3), Math.Sqrt(3) },
                { -Math.Sqrt(3), Math.Sqrt(3) }
            };

            shapeFunctionsAtNodes = new double[numNodes, numGPs];
            for (int node = 0; node < numNodes; ++node)
            {
                double r = nodeCoordinates[node, 0];
                double s = nodeCoordinates[node, 1];
                for (int gp = 0; gp < numGPs; ++gp)
                {
                    shapeFunctionsAtNodes[node, gp] = shapeFunctions[gp](r, s);
                }
            }
        }

        public IReadOnlyList<GaussPoint2D> GaussPoints { get { return gaussPoints; } }

        public Tensor2D ExtrapolateTensorFromGaussPoints(IReadOnlyList<Tensor2D> tensorsAtGPs, INaturalPoint2D point)
        {
            double r = Math.Sqrt(3) * point.Xi;
            double s = Math.Sqrt(3) * point.Eta;
            double tensorXX = 0, tensorYY = 0, tensorXY = 0;
            for (int gp = 0; gp < numGPs; ++gp)
            { 
                double shapeValue = shapeFunctions[gp](r, s);
                tensorXX += shapeValue * tensorsAtGPs[gp].XX;
                tensorYY += shapeValue * tensorsAtGPs[gp].YY;
                tensorXY += shapeValue * tensorsAtGPs[gp].XY;
            }
            return new Tensor2D(tensorXX, tensorYY, tensorXY);
        }

        public IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGPs)
        {
            var nodalTensors = new Tensor2D[numNodes];
            for (int node = 0; node < numNodes; ++node)
            {
                double tensorXX = 0, tensorYY = 0, tensorXY = 0;
                for (int gp = 0; gp < numGPs; ++gp)
                {
                    double shapeValue = shapeFunctionsAtNodes[node, gp];
                    tensorXX += shapeValue * tensorsAtGPs[gp].XX;
                    tensorYY += shapeValue * tensorsAtGPs[gp].YY;
                    tensorXY += shapeValue * tensorsAtGPs[gp].XY;
                }
                nodalTensors[node] = new Tensor2D(tensorXX, tensorYY, tensorXY);
            }
            return nodalTensors;
        }

        public Vector2 ExtrapolateVectorFromGaussPoints(IReadOnlyList<Vector2> vectorsAtGPs, INaturalPoint2D point)
        {
            double r = Math.Sqrt(3) * point.Xi;
            double s = Math.Sqrt(3) * point.Eta;
            var vector = new double[2];
            for (int gp = 0; gp < numGPs; ++gp)
            {
                double shapeValue = shapeFunctions[gp](r, s);
                vector[0] += shapeValue * vectorsAtGPs[gp][0];
                vector[1] += shapeValue * vectorsAtGPs[gp][1];
            }
            return Vector2.CreateFromArray(vector);
        }

        public IReadOnlyList<Vector2> ExtrapolateVectorFromGaussPointsToNodes(IReadOnlyList<Vector2> vectorsAtGPs)
        {
            var nodalVectors = new Vector2[numNodes];
            for (int node = 0; node < numNodes; ++node)
            {
                var vector = new double[2];
                for (int gp = 0; gp < numGPs; ++gp)
                {
                    double shapeValue = shapeFunctionsAtNodes[node, gp];
                    vector[0] += shapeValue * vectorsAtGPs[gp][0];
                    vector[1] += shapeValue * vectorsAtGPs[gp][1];
                }
                nodalVectors[node] = Vector2.CreateFromArray(vector);
            }
            return nodalVectors;
        }
    }
}
