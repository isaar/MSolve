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
    class MultipleCracksDisjoint : ICrackDescription
    {
        private readonly IReadOnlyList<TrackingExteriorCrackLSM> cracks;

        public MultipleCracksDisjoint(IReadOnlyList<TrackingExteriorCrackLSM> cracks)
        {
            this.cracks = cracks;
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesAll
        {
            get
            {
                var union = new Dictionary<CrackBodyEnrichment2D, ISet<XNode2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var enrichementNodes in crack.CrackBodyNodesAll)
                    {
                        union.Add(enrichementNodes.Key, enrichementNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesNew
        {
            get
            {
                var union = new Dictionary<CrackBodyEnrichment2D, ISet<XNode2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var enrichementNodes in crack.CrackBodyNodesNew)
                    {
                        union.Add(enrichementNodes.Key, enrichementNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesModified
        {
            get
            {
                var union = new Dictionary<CrackBodyEnrichment2D, ISet<XNode2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var enrichementNodes in crack.CrackBodyNodesModified)
                    {
                        union.Add(enrichementNodes.Key, enrichementNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesNearModified
        {
            get
            {
                var union = new Dictionary<CrackBodyEnrichment2D, ISet<XNode2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var enrichementNodes in crack.CrackBodyNodesNearModified)
                    {
                        union.Add(enrichementNodes.Key, enrichementNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesRejected
        {
            get
            {
                var union = new Dictionary<CrackBodyEnrichment2D, ISet<XNode2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var enrichementNodes in crack.CrackBodyNodesRejected)
                    {
                        union.Add(enrichementNodes.Key, enrichementNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode2D>> CrackTipNodesNew
        {
            get
            {
                var union = new Dictionary<CrackTipEnrichments2D, ISet<XNode2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var enrichementNodes in crack.CrackTipNodesNew)
                    {
                        union.Add(enrichementNodes.Key, enrichementNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode2D>> CrackTipNodesOld
        {
            get
            {
                var union = new Dictionary<CrackTipEnrichments2D, ISet<XNode2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var enrichementNodes in crack.CrackTipNodesOld)
                    {
                        union.Add(enrichementNodes.Key, enrichementNodes.Value);
                    }
                }
                return union;
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

        public IReadOnlyList<ISingleCrack> SingleCracks { get {return cracks;} }

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
