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

//TODO: replace the boilerplate code with one or more (or generic) Union() private methods
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
                    foreach (var enrichmentNodes in crack.CrackBodyNodesAll)
                    {
                        union.Add(enrichmentNodes.Key, enrichmentNodes.Value);
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
                    foreach (var enrichmentNodes in crack.CrackBodyNodesNew)
                    {
                        union.Add(enrichmentNodes.Key, enrichmentNodes.Value);
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
                    foreach (var enrichmentNodes in crack.CrackBodyNodesModified)
                    {
                        union.Add(enrichmentNodes.Key, enrichmentNodes.Value);
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
                    foreach (var enrichmentNodes in crack.CrackBodyNodesNearModified)
                    {
                        union.Add(enrichmentNodes.Key, enrichmentNodes.Value);
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
                    foreach (var enrichmentNodes in crack.CrackBodyNodesRejected)
                    {
                        union.Add(enrichmentNodes.Key, enrichmentNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyList<ICartesianPoint2D> CrackTips
        {
            get
            {
                var tips = new List<ICartesianPoint2D>();
                foreach (var crack in cracks) tips.AddRange(crack.CrackTips);
                return tips;
            }
        }

        public IReadOnlyDictionary<ICartesianPoint2D, IReadOnlyList<XContinuumElement2D>> CrackTipElements
        {
            get
            {
                var union = new Dictionary<ICartesianPoint2D, IReadOnlyList<XContinuumElement2D>>();
                foreach (var crack in cracks)
                {
                    foreach (var tipElements in crack.CrackTipElements)
                    {
                        union.Add(tipElements.Key, tipElements.Value);
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
                    foreach (var enrichmentNodes in crack.CrackTipNodesNew)
                    {
                        union.Add(enrichmentNodes.Key, enrichmentNodes.Value);
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
                    foreach (var enrichmentNodes in crack.CrackTipNodesOld)
                    {
                        union.Add(enrichmentNodes.Key, enrichmentNodes.Value);
                    }
                }
                return union;
            }
        }

        public IReadOnlyDictionary<ICartesianPoint2D, IPropagator> CrackTipPropagators
        {
            get
            {
                var union = new Dictionary<ICartesianPoint2D, IPropagator>();
                foreach (var crack in cracks)
                {
                    foreach (var tipPropagator in crack.CrackTipPropagators)
                    {
                        union.Add(tipPropagator.Key, tipPropagator.Value);
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
