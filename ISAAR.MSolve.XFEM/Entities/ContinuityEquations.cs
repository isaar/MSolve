using ISAAR.MSolve.XFEM.FreedomDegrees;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.XFEM.Entities
{
    class ContinuityEquations
    {
        public ContinuityEquations(int numEquations, IReadOnlyList<XNode> nodes, IReadOnlyList<EnrichedDof> dofs,
            IReadOnlyList<XSubdomain2D_old> positiveSubdomains, IReadOnlyList<XSubdomain2D_old> negativeSubdomains)
        {
            this.NumEquations = numEquations;
            this.Nodes = nodes;
            this.Dofs = dofs;
            this.PositiveSubdomains = positiveSubdomains;
            this.NegativeSubdomains = negativeSubdomains;
        }

        public enum EquationOrder { DofMajor, SubdomainMajor};

        public int NumEquations { get; }
        public IReadOnlyList<XNode> Nodes { get; }
        public IReadOnlyList<EnrichedDof> Dofs { get; }
        public IReadOnlyList<XSubdomain2D_old> PositiveSubdomains { get; }
        public IReadOnlyList<XSubdomain2D_old> NegativeSubdomains { get; }


        public class Builder
        {
            private readonly Dictionary<XNode, Dictionary<EnrichedDof, SortedSet<XSubdomain2D_old>>> data;

            public Builder()
            {
                this.data = new Dictionary<XNode, Dictionary<EnrichedDof, SortedSet<XSubdomain2D_old>>>();
            }

            public ContinuityEquations Build(EquationOrder order)
            {
                var nodes = new List<XNode>();
                var dofs = new List<EnrichedDof>();
                var positiveSubdomains = new List<XSubdomain2D_old>();
                var negativeSubdomains = new List<XSubdomain2D_old>();

                if (order == EquationOrder.DofMajor)
                {
                    foreach (var nodeData in data)
                    {
                        XNode node = nodeData.Key;
                        foreach (var dofData in nodeData.Value)
                        {
                            EnrichedDof dof = dofData.Key;
                            if (dofData.Value.Count > 1) // There must be >= 2 subdomains, in order to enforce continuity
                            {
                                XSubdomain2D_old[] subdomains = dofData.Value.ToArray();
                                for (int i = 0; i < subdomains.Length - 1; ++i)
                                {
                                    nodes.Add(node);
                                    dofs.Add(dof);
                                    positiveSubdomains.Add(subdomains[i]);
                                    negativeSubdomains.Add(subdomains[i + 1]);
                                }
                            }
                        }
                    }
                    return new ContinuityEquations(nodes.Count, nodes, dofs, positiveSubdomains, negativeSubdomains);
                }
                else
                {
                    throw new NotImplementedException();
                }
            }

            public void Register(XNode node, EnrichedDof dof, XSubdomain2D_old subdomain)
            {
                bool nodeExists = data.TryGetValue(node, out Dictionary<EnrichedDof, SortedSet<XSubdomain2D_old>> nodeData);
                if (!nodeExists)
                {
                    nodeData = new Dictionary<EnrichedDof, SortedSet<XSubdomain2D_old>>();
                    data.Add(node, nodeData);
                }

                bool dofExists = nodeData.TryGetValue(dof, out SortedSet<XSubdomain2D_old> subdomains);
                if (!dofExists)
                {
                    subdomains = new SortedSet<XSubdomain2D_old>();
                    nodeData.Add(dof, subdomains);
                }

                subdomains.Add(subdomain);
            }

        }
    }
}
