using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Tensors;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Interpolation.GaussPointSystems
{
    interface IGaussPointSystem
    {
        IReadOnlyList<GaussPoint> GaussPoints { get; }

        Tensor2D ExtrapolateTensorFromGaussPoints(IReadOnlyList<Tensor2D> tensorsAtGPs, NaturalPoint point);
        IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGPs);
        Vector2 ExtrapolateVectorFromGaussPoints(IReadOnlyList<Vector2> vectorsAtGPs, NaturalPoint point);
        IReadOnlyList<Vector2> ExtrapolateVectorFromGaussPointsToNodes(IReadOnlyList<Vector2> vectorsAtGPs);
    }
}