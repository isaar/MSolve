using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    interface IBenchmark
    {
        string Name { get; }

        ICrackDescription Crack { get; }
        IDomainDecomposer Decomposer { get; }

        Model2D Model { get; }
        string PlotDirectory { get; }
        Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> PossibleEnrichments { get; }

        void Analyze(ISolver solver);
        void InitializeModel();
    }
}
