using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Tensors;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Interpolation.GaussPointSystems
{
    interface IGaussPointSystem
    {
        IReadOnlyList<GaussPoint2D> GaussPoints { get; }

        Tensor2D ExtrapolateTensorFromGaussPoints(IReadOnlyList<Tensor2D> tensorsAtGPs, NaturalPoint2D point);
        IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGPs);
        Vector2 ExtrapolateVectorFromGaussPoints(IReadOnlyList<Vector2> vectorsAtGPs, NaturalPoint2D point);
        IReadOnlyList<Vector2> ExtrapolateVectorFromGaussPointsToNodes(IReadOnlyList<Vector2> vectorsAtGPs);
    }
}