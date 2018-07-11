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
    interface IGaussPointSystem
    {
        IReadOnlyList<GaussPoint2D> GaussPoints { get; }

        Tensor2D ExtrapolateTensorFromGaussPoints(IReadOnlyList<Tensor2D> tensorsAtGPs, INaturalPoint2D point);
        IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGPs);
        Vector2 ExtrapolateVectorFromGaussPoints(IReadOnlyList<Vector2> vectorsAtGPs, INaturalPoint2D point);
        IReadOnlyList<Vector2> ExtrapolateVectorFromGaussPointsToNodes(IReadOnlyList<Vector2> vectorsAtGPs);
    }
}