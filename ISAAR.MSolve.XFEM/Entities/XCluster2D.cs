using System;
using System.Collections.Generic;
using System.Text;

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

        public void AddSubdomain(XSubdomain2D subdomain)
        {
            subdomains.Add(subdomain);
        }
    }
}
