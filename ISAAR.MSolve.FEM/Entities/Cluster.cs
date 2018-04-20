using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class Cluster
    {
        private readonly IList<Subdomain> subdomains = new List<Subdomain>();

        public IList<Subdomain> Subdomains
        {
            get { return subdomains; }
        }

        public int ID { get; set; }
    }
}
