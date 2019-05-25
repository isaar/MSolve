using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019
{
    public static class Utilities
    {
        public static Dictionary<int, INode[]> FindCornerNodesFromCrosspoints2D(IStructuralModel model)
        {
            //TODO: This is also done by the analyzer. Perhaps it should not.
            model.ConnectDataStructures();
            var cornerNodes = new Dictionary<int, INode[]>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                // In 2D, if multiplicity > 2, the node is a crosspoint 
                cornerNodes[subdomain.ID] = subdomain.Nodes.Where(node => node.SubdomainsDictionary.Count > 2).ToArray();
            }
            return cornerNodes;
        }
    }
}
