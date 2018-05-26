using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    interface ICrackDescription
    {
        ISet<XNode2D> CrackBodyNodesAll { get; }
        ISet<XNode2D> CrackBodyNodesNew { get; }
        ISet<XNode2D> CrackBodyNodesModified { get; }
        ISet<XNode2D> CrackBodyNodesNearModified { get;}
        ISet<XNode2D> CrackTipNodesNew { get; }
        ISet<XNode2D> CrackTipNodesOld { get; }
        IReadOnlyList<EnrichedDof> DofsHeaviside { get; }
        IReadOnlyList<EnrichedDof> DofsTip { get; }
        ISet<XContinuumElement2D> ElementsModified { get; }
        IReadOnlyList<IEnrichmentItem2D> Enrichments { get; }
        BiMesh2D Mesh { get; } //TODO: abstract this

        IReadOnlyList<ICartesianPoint2D> GetCrackTips();
        IReadOnlyList<IPropagator> GetCrackTipPropagators();
        void Propagate(IDofOrderer dofOrderer, Vector totalFreeDisplacements, Vector totalConstrainedDisplacements);
        void UpdateEnrichments();
    }
}
