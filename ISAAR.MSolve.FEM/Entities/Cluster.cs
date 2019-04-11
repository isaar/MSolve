using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class Cluster
    {
        private readonly IList<Subdomain_v2> subdomains = new List<Subdomain_v2>();

        public IList<Subdomain_v2> Subdomains
        {
            get { return subdomains; }
        }

        public int ID { get; set; }
    }
}
