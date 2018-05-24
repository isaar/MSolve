using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    interface IBenchmark
    {
        TrackingExteriorCrackLSM Crack { get; }
        IDecomposer Decomposer { get; }
        IReadOnlyList<XNode2D> EnrichedArea { get; }
        Model2D Model { get; }
        IReadOnlyList<ICartesianPoint2D> Analyze(ISolver solver);
        void InitializeModel();
    }
}
