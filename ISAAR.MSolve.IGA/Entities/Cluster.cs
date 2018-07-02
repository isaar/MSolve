using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Entities
{
    public class Cluster
    {
        private readonly IList<Patch> patches = new List<Patch>();

        public IList<Patch> Patches
        {
            get { return patches; }
        }

        public int ID { get; set; }
    }

    
}
