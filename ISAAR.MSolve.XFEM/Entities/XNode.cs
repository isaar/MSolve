using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    public class XNode : Node
    {
        public XNode(int id, double x, double y) : base(id, x, y)
        {
            this.EnrichmentItems = new Dictionary<IEnrichmentItem2D, double[]>();
        }

        public XNode(int id, double x, double y, double z) : base(id, x, y, z)
        {
            this.EnrichmentItems = new Dictionary<IEnrichmentItem2D, double[]>();
        }

        public Dictionary<IEnrichmentItem2D, double[]> EnrichmentItems { get; }
        public bool IsEnriched { get { return EnrichmentItems.Count > 0; } }

        public int EnrichedDofsCount
        {
            get
            {
                int count = 0;
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys) count += enrichment.Dofs.Count;
                return count;
            }
        }

        public IReadOnlyList<EnrichedDof> EnrichedDofs
        {
            get
            {
                var dofs = new List<EnrichedDof>();
                foreach (IEnrichmentItem2D enrichment in EnrichmentItems.Keys) dofs.AddRange(enrichment.Dofs);
                return dofs;
            }
        }

        //TODO: Redesign these and their counterparts in Node. Connectivity should be done using the Dircetization interfaces.
        #region connectivity
        public new Dictionary<int, IXFiniteElement> ElementsDictionary { get; } = new Dictionary<int, IXFiniteElement>();
        //public new Dictionary<int, XSubdomain> SubdomainsDictionary { get; } = new Dictionary<int, XSubdomain>();

        public void BuildXSubdomainDictionary()
        {
            foreach (IXFiniteElement element in ElementsDictionary.Values)
            {
                SubdomainsDictionary[element.Subdomain.ID] = element.Subdomain;
            }
        }
        #endregion
    }
}