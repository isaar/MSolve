using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes
{
    public class UsedDefinedCornerNodes : ICornerNodeSelection
    {
        private readonly Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains;

        public UsedDefinedCornerNodes(Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains)
        {
            this.cornerNodesOfSubdomains = cornerNodesOfSubdomains;
        }

        public Dictionary<int, HashSet<INode>> SelectCornerNodesOfSubdomains()
            => cornerNodesOfSubdomains;
    }
}
