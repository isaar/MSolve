using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// Inter-subdomain communication 
    /// </summary>
    class XCluster2D
    {
        private readonly List<XSubdomain2D> subdomains;

        public XCluster2D()
        {
            this.subdomains = new List<XSubdomain2D>();
        }

        public XCluster2D(IReadOnlyList<XSubdomain2D> subdomains)
        {
            this.subdomains = new List<XSubdomain2D>(subdomains);
        }

        public XClusterDofOrderer DofOrderer { get; private set; }
        public IReadOnlyList<XSubdomain2D> Subdomains { get { return subdomains; } }

        public void AddSubdomain(XSubdomain2D subdomain)
        {
            subdomains.Add(subdomain);
        }

        public XSubdomain2D FindSubdomainOfElement(XContinuumElement2D element)
        {
            foreach (var subdomain in subdomains)
            {
                if (subdomain.Elements.Contains(element)) return subdomain;
            }
            throw new KeyNotFoundException("This element does not belong to any subdomain.");
        }

        public void OrderDofs(Model2D model)
        {
            DofOrderer = XClusterDofOrderer.CreateNodeMajor(model, this);
        }
    }
}
