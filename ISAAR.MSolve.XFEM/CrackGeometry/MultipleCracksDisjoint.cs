using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Interpolation;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    class MultipleCracksDisjoint: ICrackDescription
    {
        private readonly IReadOnlyList<TrackingExteriorCrackLSM> cracks;

        public MultipleCracksDisjoint(IReadOnlyList<TrackingExteriorCrackLSM> cracks)
        {
            this.cracks = cracks;
        }

        public ISet<XNode2D> CrackBodyNodesAll
        {
            get
            {
                var nodes = new HashSet<XNode2D>();
                foreach (var crack in cracks) nodes.UnionWith(crack.CrackBodyNodesAll);
                return nodes;
            }
        }

        public ISet<XNode2D> CrackBodyNodesNew
        {
            get
            {
                var nodes = new HashSet<XNode2D>();
                foreach (var crack in cracks) nodes.UnionWith(crack.CrackBodyNodesNew);
                return nodes;
            }
        }

        public ISet<XNode2D> CrackBodyNodesModified
        {
            get
            {
                var nodes = new HashSet<XNode2D>();
                foreach (var crack in cracks) nodes.UnionWith(crack.CrackBodyNodesModified);
                return nodes;
            }
        }

        public ISet<XNode2D> CrackBodyNodesNearModified
        {
            get
            {
                var nodes = new HashSet<XNode2D>();
                foreach (var crack in cracks) nodes.UnionWith(crack.CrackBodyNodesNearModified);
                return nodes;
            }
        }

        public ISet<XNode2D> CrackTipNodesNew
        {
            get
            {
                var nodes = new HashSet<XNode2D>();
                foreach (var crack in cracks) nodes.UnionWith(crack.CrackTipNodesNew);
                return nodes;
            }
        }

        public ISet<XNode2D> CrackTipNodesOld
        {
            get
            {
                var nodes = new HashSet<XNode2D>();
                foreach (var crack in cracks) nodes.UnionWith(crack.CrackTipNodesOld);
                return nodes;
            }
        }

        public IReadOnlyList<EnrichedDof> DofsHeaviside
        {
            get
            {
                var dofs = new List<EnrichedDof>();
                foreach (var crack in cracks) dofs.AddRange(crack.DofsHeaviside);
                return dofs;
            }
        }

        public IReadOnlyList<EnrichedDof> DofsTip
        {
            get
            {
                var dofs = new List<EnrichedDof>();
                foreach (var crack in cracks) dofs.AddRange(crack.DofsTip);
                return dofs;
            }
        }

        public ISet<XContinuumElement2D> ElementsModified
        {
            get
            {
                var elements = new HashSet<XContinuumElement2D>();
                foreach (var crack in cracks) elements.UnionWith(crack.ElementsModified);
                return elements;
            }
        }

        public BiMesh2D Mesh { get { return cracks[0].Mesh; } } //TODO: during construction check that all cracks use the same mesh

        public IReadOnlyList<IEnrichmentItem2D> Enrichments
        {
            get
            {
                var enrichments = new List<IEnrichmentItem2D>();
                foreach (var crack in cracks) enrichments.AddRange(crack.Enrichments);
                return enrichments;
            }
        }

        public IReadOnlyList<IPropagator> GetCrackTipPropagators()
        {
            var propagators = new List<IPropagator>();
            foreach (var crack in cracks) propagators.AddRange(crack.GetCrackTipPropagators());
            return propagators;
        }

        public IReadOnlyList<ICartesianPoint2D> GetCrackTips()
        {
            var tips = new List<ICartesianPoint2D>();
            foreach (var crack in cracks) tips.AddRange(crack.GetCrackTips());
            return tips;
        }

        public void Propagate(IDofOrderer dofOrderer, Vector totalFreeDisplacements, Vector totalConstrainedDisplacements)
        {
            foreach (var crack in cracks) crack.Propagate(dofOrderer, totalFreeDisplacements, totalConstrainedDisplacements);
        }

        public void UpdateEnrichments()
        {
            foreach (var crack in cracks) crack.UpdateEnrichments();
        }
    }
}
