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
    class Tri3GPSystem: IGaussPointSystem
    {
        private const int numGPs = 3;
        private const int numNodes = 3; //TODO: Actually this is used for Tri6

        // The order is important
        private static readonly GaussPoint2D[] gaussPoints = new GaussPoint2D[]
        {
            new GaussPoint2D(1.0/6.0, 1.0/6.0, 1.0/6.0),
            new GaussPoint2D(2.0/3.0, 1.0/6.0, 1.0/6.0),
            new GaussPoint2D(1.0/6.0, 2.0/3.0, 1.0/6.0)
        };

        private static readonly double[,] shapeFunctionsAtNodes;

        private static readonly Func<double, double, double>[] shapeFunctions = new Func<double, double, double>[]
        {
            (r, s) => { return 1 - r - s; },
            (r, s) => { return r; },
            (r, s) => { return s; }
        };

        static Tri3GPSystem()
        {
            // Coordinates of nodes in the GP system: r=2*xi-1/3, s=2*eta-1/3
            double[,] nodeCoordinates = new double[,]
            {
                { -1.0/3.0, -1.0/3.0 },
                { 5.0/3.0, -1.0/3.0 },
                { -1.0/3.0, 5.0/3.0 }
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
            double r = 2.0 * point.Xi - 1.0 / 3.0;
            double s = 2.0 * point.Eta - 1.0 / 3.0;
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
            double r = 2.0 * point.Xi - 1.0 / 3.0;
            double s = 2.0 * point.Eta - 1.0 / 3.0;
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
