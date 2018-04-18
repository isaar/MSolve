using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Interpolation.GaussPointSystems
{
    interface IGaussPointSystem
    {
        IReadOnlyList<INaturalPoint2D> GaussPoints { get; }

        IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGPs);
    }
}