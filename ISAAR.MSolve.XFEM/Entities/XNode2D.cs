using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    public class XNode2D : CartesianPoint2D, INode
    {
        private readonly int id;

        public XNode2D(int id, double x, double y) : base(x, y)
        {
            this.id = id;
            this.EnrichmentItems = new Dictionary<IEnrichmentItem2D, double[]>();
        }

        public Dictionary<IEnrichmentItem2D, double[]> EnrichmentItems { get; }
        public bool IsEnriched { get { return EnrichmentItems.Count > 0; } }

        public int EnrichedDofsCount
        {
            get
            {
                int count = 0;
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys)
                {
                    count += enrichment.Dofs.Count;
                }
                return count;
            }
        }

        public ISet<EnrichedDof> EnrichedDofs
        {
            get
            {
                var dofs = new HashSet<EnrichedDof>();
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys)
                {
                    dofs.UnionWith(enrichment.Dofs);
                }
                return dofs;
            }
        }

        public int ID { get => id; set => throw new NotImplementedException("Cannot be modified"); }
        public double Z { get => throw new NotImplementedException("2D node"); set => throw new NotImplementedException("2D node"); }

        public List<Constraint> Constraints => throw new NotImplementedException();

        public Dictionary<int, ISubdomain> SubdomainsDictionary => throw new NotImplementedException();

        double INode.X { get => base.X; set => throw new NotImplementedException("Cannot be modified"); }
        double INode.Y { get => base.Y; set => throw new NotImplementedException("Cannot be modified"); }

        public int CompareTo(INode other) => this.ID - other.ID;
    }
}